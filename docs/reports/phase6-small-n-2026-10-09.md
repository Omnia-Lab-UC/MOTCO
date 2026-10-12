# Phase 6 small-n operating study at the SEA-AD design — 2026-10-09

> Follows `examples/trajectory_power_study/phase5_report_template.md` section for section. The study is split
> into two profiles (spec "Phase 6 small-n study profiles", split rule), so each section reports both runs.
> Rates and Monte Carlo SEs are read from the two runs' `report/` directories. The few figures that are not in
> `report/` (selected-rank breakdowns, anchor reproduction) are computed from `merged.jsonl`, whose sha256 is
> in each `PROVENANCE.txt`, and are marked as such.

**Runs:** `results/phase6-small-n-magnitude-2026-10-09/` (magnitude) and `results/phase6-small-n-2026-10-09/`
(orientation, shape, translation) · **Configs:** `examples/trajectory_power_study/phase6_small_n_magnitude.json`
(sha256 `41211d61…`) and `examples/trajectory_power_study/phase6_small_n_study.json` (sha256 `18cf0fc7…`) ·
**Gate decisions:** magnitude **PROCEED**, orientation/shape **HOLD** (advisory, design D6; from each run's
`report/phase4_gate_decision.json`).

**In one paragraph.** At the SEA-AD MTG astrocyte design (n = 80, three stages, cells of 9–28, two blocks),
the three RRPP tests hold their size. Every Type I cell of the baseline column is at or below α (0.010–0.040, bound 0.0695), and the three-block columns stay inside the bound. The tests
differ sharply in power. `delta` detects a joint size change of 10 % about half the time (0.542) and of 25 %
almost always (0.956), and the magnitude construction leaves `angle` and `shape` at the floor. `angle` and
`shape` stay weak even at the strongest orientation and shape constructions the generator can build: 0.300 and
0.228 at e = 1.00, where every differential site has been relocated. They are also not specific. The orientation
construction rejects `shape` (0.464) and `delta` (0.354) more often than `angle`, and the shape construction
rejects `delta` (0.330) more often than `shape`. For the case study, a `delta` rejection is interpretable
evidence that the trajectories differ. A `delta` non-rejection is evidence against size changes of about 25 % or
more. `angle` and `shape` rejections are interpretable as evidence of a difference, but not of which kind, and
their non-rejections carry no information at n = 80 (section 8).

## 1. Configuration and provenance

**Design point.** The study reproduces the SEA-AD MTG astrocyte case-study design fixed on 2026-10-09. Sex is
the group (A = F, B = M; the surgery is applied to the smaller group, the conservative assignment, design D4).
ADNC is merged into three stages (Not AD + Low / Intermediate / High). The donor table at ≥30 nuclei per donor
per omic is `group_stage_sizes = [[11, 10, 28], [9, 10, 12]]` (n = 80). All three InterSIM blocks are generated,
but only methylation + expression are measured (`integration_params.layers`), as the two-block analogue of
ATAC + RNA (design D2). Every trajectory surgery acts on methylation indicators, so each mode's construction is
present in the measured space. Measurement is pooled PLS on M-value methylation with the committed
stage-supervised double-CV rank, and `magnitude_kind = joint`. `p_dmp = 0.1` and ρ = 0 are copied from
Phase 5. Only the stage count and sizes change.

**Design grid.** Both profiles cross the ≥30-nuclei and ≥50-nuclei tables (`[[10, 9, 25], [9, 9, 12]]`, n = 74)
with two-block and three-block measurement (`layers` null), giving four columns. The cohort-size contrast is the
baseline vs the ≥50-nuclei column at two blocks. The block-count contrast is the baseline vs the three-block
column at n = 80. Columns that differ only in `layers` evaluate the same generated datasets. The fourth column
(≥50 nuclei, three blocks) exists because the engine crosses design-grid axes (design D4, accepted 2026-10-09).

**Effect axes.** These come from the pilot (`results/phase6-small-n-pilot-2026-10-09/NOTES.md`; each profile's
`metadata.effect_axis` gives the reason for each point). Magnitude: 0 / 0.02 / 0.05 / 0.10 / 0.25 / 0.50 / 1.00,
covering the measured rise and saturation. Orientation, shape and translation: 0 / 0.25 / 0.50 / 0.75 / 1.00.
For orientation and shape, e = 1.00 is the construction maximum, because `_relocate_rows` clamps the relocated
fraction at 1.

**Monte Carlo sizing and seeds.** 500 replicates × 999 permutations per cell. Base seed 800 and the matched-seed
family `phase6-small-n` are shared with the pilot. This makes the shared zero-effect anchor the same datasets in
the pilot and both profiles, and the family is independent of every Phase 5 record. Both configs derive from
`phase6_small_n_pilot.json` (`metadata.derives_from`) and name each other as split siblings.

**Provenance.** See each run's `PROVENANCE.txt`; the fields the template requires are all present. In summary:
- Records were produced at code revision `5325f6d`; the reports were rendered at `946e471`, which adds only a
  reporting fix (section 2).
- The runs used cluster `ing`, partition `512x1024`, AMD EPYC 7662, with `--exclude=n3` and later `n3,n4`.
  `STUDY_N_JOBS` was unset and BLAS was pinned to one thread.
- Software: python 3.11.16, motco 0.6.0, numpy 2.3.5, scikit-learn 1.8.0, scipy 1.16.3, uv 0.12.10.
- Layout: 100 shards per profile, error policy `record`, and `n_jobs = 1` in every record (each
  `report/report_contract.json`).
- Magnitude: 30 cells × 500 = 15,000 units, median 11.6 s per unit, 47.6 core-hours.
- Orientation/shape/translation: 54 cells × 500 = 27,000 units, median 13.2 s per unit, 102.4 core-hours.

## 2. Unit and failure accounting

| | magnitude | orientation / shape / translation |
|---|---|---|
| units expected | 15,000 (30 cells × 500) | 27,000 (54 cells × 500) |
| units present exactly once | 15,000 | 27,000 |
| cells present with exactly 500 replicates | 30 / 30 | 54 / 54 |
| parameter-signature mismatches | 0 | 0 |
| `status = failed` records | 0 | 0 |
| censored surgeries (`report/realized_surgery.csv`, `censored_fraction`) | 0 at every cell | 0 at every cell |
| attribution (`report/driver_report.csv`) | not requested (15,000) | 4 orientation cells: 500 eligible / 500 computed / 0 failed each; 25,000 not requested |

**Node incident.** At 19:40Z node n4 ran out of memory under other users' jobs, and the 12 magnitude shards
running there stalled without writing. They were cancelled and requeued, and no unit was lost or duplicated
(magnitude `PROVENANCE.txt`, `node_incident`).

**Report fix.** The first `report` at `5325f6d` crashed before writing anything: design-point coordinates for a
sequence-valued axis come back from JSONL as unhashable lists. `946e471` fixes this and adds a regression test
(`tests/test_study_design_point_report.py::test_list_valued_design_axes_survive_a_jsonl_round_trip`). Records
and statistics are unaffected.

## 3. Gate decision

The Phase 4 gate is declared as **advisory** acceptance (design D6). At n = 80 a missed 0.80 power floor is a
finding about the case study's design, not a defect.

**Magnitude profile — PROCEED** (`report/phase4_gate_decision.json`; rules reduced to magnitude as in the Phase 5
re-measurement):
- Mandatory power: magnitude/`delta` reaches 1.000 at e = 1.00 against the 0.80 floor.
- Mandatory controls at e = 1.00, bound α + 2·SE = 0.0695: magnitude/`angle` 0.008, magnitude/`shape` 0.002.
- Type I inflation: anchor 0.026 / 0.024 / 0.020 (`delta` / `angle` / `shape`); `type_i_baseline` `none`
  0.016 / 0.028 / 0.018.
- Completeness met.

**Orientation/shape profile — HOLD** (`report/phase4_gate_decision.json`). Two mandatory power rules fail, and
both curves are monotone (`tolerated_monotone=True`):
- orientation/`angle`: top rate 0.300 against the 0.80 floor.
- shape/`shape`: top rate 0.228 against the 0.80 floor.

Every Type I inflation observation passes (`report/phase4_gate.csv`): the anchor, `type_i_baseline` `none` and
translation, and translation at every e (0.010–0.040, bound 0.0695). The descriptive roles (orientation/`delta`,
orientation/`shape`, shape/`delta`, shape/`angle`) are reported in section 6, and completeness is met. No
confirmation re-run is required: a HOLD on mandatory power at the top effect is not a marginal exceedance. The
HOLD is the measured answer to "can orientation and shape be detected at n = 80?".

## 4. Type I error

Null claims read the `type_i_baseline` cells, a separate seed family from the power grid (`report/type_i_table.csv`):

| cell | `delta` | `angle` | `shape` |
|---|---|---|---|
| `none` (`type_i_baseline-ca766c21680e`) | 0.016 ± 0.006 | 0.028 ± 0.007 | 0.018 ± 0.006 |
| translation at e = 1.00 (`type_i_baseline-74c0457435ca`) | 0.034 ± 0.008 | 0.040 ± 0.009 | 0.018 ± 0.006 |

Both cells have the same cell ids and signatures in the two profiles and are the same measurement. They are
counted once.

The shared zero-effect anchor reads 0.026 ± 0.007 / 0.024 ± 0.007 / 0.020 ± 0.006. It is
`power_primary-6f04726ce01b` in the magnitude profile and `power_primary-1a2e4088b485` in the
orientation/shape profile (`report/report_contract.json` resolves it to `[magnitude]` and to
`[orientation, shape, translation]`, `counted_as: 1`). The two cells hold the same 500 datasets. Generator seeds
and selected ranks are identical at all 500 replicate indices, observed statistics agree to within 1.4e-11, and
p-values are identical (computed from `merged.jsonl`). The anchor is therefore **one** measurement, not one per
mode or per profile. Its first 100 replicates also reproduce the pilot's anchor to within 7.3e-12 on the
observed statistics, with identical seeds and ranks (p-values differ: 199 vs 999 permutations).

Translation at e > 0, read as a negative control on all three statistics (`report/power_curves.csv`):
0.014–0.032 for `delta`, 0.018–0.030 for `angle` and 0.010–0.030 for `shape`, at every e.

The two-block measurement is **conservative** at n = 80. The baseline column's null rates sit at half of α. At
three blocks the anchor rates rise to 0.062 ± 0.011 (`delta`) and 0.056 ± 0.010 (`shape`), still inside the
bound; at ≥50 nuclei and two blocks they fall to 0.010 and 0.008 (`report/design_point_operating.csv`).

## 5. Power per mode

Baseline column (≥30 nuclei, two blocks), `report/power_curves.csv`, rate ± Monte Carlo SE:

| mode | e | `delta` | `angle` | `shape` |
|---|---|---|---|---|
| anchor | 0 | 0.026 ± 0.007 | 0.024 ± 0.007 | 0.020 ± 0.006 |
| magnitude | 0.02 | **0.052 ± 0.010** | 0.024 ± 0.007 | 0.022 ± 0.007 |
| | 0.05 | **0.206 ± 0.018** | 0.020 ± 0.006 | 0.024 ± 0.007 |
| | 0.10 | **0.542 ± 0.022** | 0.018 ± 0.006 | 0.022 ± 0.007 |
| | 0.25 | **0.956 ± 0.009** | 0.016 ± 0.006 | 0.014 ± 0.005 |
| | 0.50 | **1.000** | 0.010 ± 0.004 | 0.006 ± 0.003 |
| | 1.00 | **1.000** | 0.008 ± 0.004 | 0.002 ± 0.002 |
| orientation | 0.25 | 0.070 ± 0.011 | **0.034 ± 0.008** | 0.068 ± 0.011 |
| | 0.50 | 0.160 ± 0.016 | **0.058 ± 0.010** | 0.196 ± 0.018 |
| | 0.75 | 0.194 ± 0.018 | **0.090 ± 0.013** | 0.320 ± 0.021 |
| | 1.00 | 0.354 ± 0.021 | **0.300 ± 0.020** | 0.464 ± 0.022 |
| shape | 0.25 | 0.104 ± 0.014 | 0.042 ± 0.009 | **0.060 ± 0.011** |
| | 0.50 | 0.216 ± 0.018 | 0.076 ± 0.012 | **0.090 ± 0.013** |
| | 0.75 | 0.292 ± 0.020 | 0.112 ± 0.014 | **0.172 ± 0.017** |
| | 1.00 | 0.330 ± 0.021 | 0.138 ± 0.015 | **0.228 ± 0.019** |

**Magnitude.** The `delta` rise sits about 10× higher on the native joint axis than at n = 1200 (Phase 5
re-measurement: 0.836 at e = 0.01, 0.996 at 0.02). e is the per-omic size ratio minus one, so `delta` has about
even odds at a 10 % joint size change and is near certain at 25 %. The pilot (100 × 199) measured 0.58 and 0.93
at the same points, agreeing within pooled SE.

**Orientation, with its recorded geometry.** The anchor's pooled relative eigengap has median 0.111 and terciles
0.082 / 0.144 (`report/design_point_operating.csv`, baseline rows). Three stage means always span at most a
2-D subspace. The small eigengap, a near-isotropic triangle whose PC1 is poorly determined, comes from the
generator's isotropic endpoint ρ = 0: each stage draws its differential sites independently, so the stage means
do not trend. It is not a property of three stages. A real cohort whose stages trend would have a larger
eigengap. The Phase 5 continuity sweep found that this does not translate into orientation power under the
current construction: the orientation surgery realizes a smaller latent contrast on a trending baseline
([design-point report](phase5-design-point-pilot-2026-09-08.md), "Continuity does not help"). The `angle` null is
correspondingly wide and variable: per-replicate q95 median 84° and IQR 115° at the anchor. At the orientation
cells:

| e | median eigengap (terciles) | `angle` null q95 median (IQR) | orientation/`angle` power |
|---|---|---|---|
| 0.25 | 0.135 (0.101 / 0.173) | 59° (109°) | 0.034 ± 0.008 |
| 0.50 | 0.160 (0.120 / 0.210) | 56° (105°) | 0.058 ± 0.010 |
| 0.75 | 0.176 (0.135 / 0.222) | 65° (114°) | 0.090 ± 0.013 |
| 1.00 | 0.167 (0.127 / 0.213) | 63° (104°) | 0.300 ± 0.020 |

Orientation `angle` power within eigengap terciles at e = 1.00 (`report/eigengap_stratified_power.csv`,
`power_primary` rows):

| tercile | eigengap range | power |
|---|---|---|
| low | 0.013–0.127 | 0.114 ± 0.025 |
| middle | 0.127–0.213 | 0.449 ± 0.038 |
| high | 0.213–0.441 | 0.337 ± 0.037 |

At e = 0.75 the same terciles give 0.060 / 0.138 / 0.072. In the low tercile, `angle` power is at most 0.114 at
every e (0.018 / 0.048 / 0.060 / 0.114).

**Rank selection.** Power at e = 1.00 also depends on the CV-selected rank (`merged.jsonl`). Elsewhere the rank is
2 or 3 in essentially every replicate. At e = 1.00 in the baseline column, 257 of 500 replicates select more than
three components, because relocating every differential site creates group-specific variation that CV
retains. Those replicates reject `angle` at 0.506, against 0.082 for the 243 at rank ≤ 3.

**Comparison with the pilots.** The Phase 5 design-point and ladder pilots measured orientation `angle` at
0.88 and 0.85 (n = 1200, four stages). The Phase 6 pilot (100 × 199) measured 0.24 at e = 1.00 and the
paper-grade run 0.300 ± 0.020, a difference of 1.3 pooled SE (pilot SE about 0.043).

**Shape.** shape/`shape` rises to 0.228 at the construction maximum, against a pilot value of 0.18. The `shape`
null is not the limiting factor (q95 near 0.1 in the pilot records). The construction is small: one interior
vertex of a three-vertex path.

### Cohort-size and block-count contrasts

`report/design_point_operating.csv`, target statistic per mode, at the points where the curves move:

| mode / statistic, e | baseline (n = 80, 2 blocks) | ≥50 nuclei (n = 74, 2 blocks) | three blocks (n = 80) | n = 74, 3 blocks |
|---|---|---|---|---|
| anchor `delta` / `angle` / `shape` | 0.026 / 0.024 / 0.020 | 0.010 / 0.022 / 0.008 | 0.062 / 0.020 / 0.056 | 0.038 / 0.024 / 0.044 |
| magnitude/`delta`, 0.05 | 0.206 ± 0.018 | 0.198 ± 0.018 | 0.438 ± 0.022 | 0.436 ± 0.022 |
| magnitude/`delta`, 0.10 | 0.542 ± 0.022 | 0.502 ± 0.022 | 0.774 ± 0.019 | 0.770 ± 0.019 |
| orientation/`angle`, 0.50 | 0.058 ± 0.010 | 0.066 ± 0.011 | 0.162 ± 0.016 | 0.144 ± 0.016 |
| orientation/`angle`, 1.00 | 0.300 ± 0.020 | 0.112 ± 0.014 | 0.294 ± 0.020 | 0.210 ± 0.018 |
| shape/`shape`, 0.50 | 0.090 ± 0.013 | 0.088 ± 0.013 | 0.306 ± 0.021 | 0.276 ± 0.020 |
| shape/`shape`, 1.00 | 0.228 ± 0.019 | 0.138 ± 0.015 | 0.472 ± 0.022 | 0.448 ± 0.022 |

**Cohort size (baseline vs ≥50 nuclei, two blocks).** Dropping six donors leaves magnitude essentially unchanged
(0.206 → 0.198; 0.542 → 0.502, within 1.3 pooled SE). It costs orientation/`angle` at e = 1.00 most of its power
(0.300 → 0.112). That drop is a rank-selection effect, not a smooth sample-size effect. At n = 74 only 39 of 500
replicates select more than three components, against 257 at n = 80, and within each rank class power is
similar (0.436 vs 0.506 above rank 3; 0.085 vs 0.082 at or below it; `merged.jsonl`). shape/`shape` at e = 1.00
falls from 0.228 to 0.138. Below the construction maximum, orientation and shape are unchanged within SE.

**Block count (baseline vs three blocks, n = 80, same datasets).** The third standardized block raises every
mode's power where the curve is rising:
- magnitude/`delta` at 0.05: 0.206 → 0.438;
- orientation/`angle` at 0.50: 0.058 → 0.162;
- shape/`shape` at 0.50: 0.090 → 0.306.

It also raises the null rates of `delta` and `shape` from about half of α to near α (0.026 → 0.062,
0.020 → 0.056), so part of the gain is the two-block measurement's conservatism being removed. Orientation at
e = 1.00 is the exception (0.300 vs 0.294): there both columns are dominated by the rank-selection effect above.

**Interaction column (n = 74, three blocks).** It follows the three-block column within SE, except
orientation/`angle` at e = 1.00 (0.210), where rank selection again carries the difference.

## 6. Cross-talk

The baseline column's off-diagonals are shown below (`report/specificity_matrix.csv` at the top effect,
`report/power_curves.csv` per effect, and `report/phase4_localization.csv` for where each first becomes
material).

**Magnitude's off-diagonals are mandatory controls.** magnitude/`angle` is 0.008–0.024 and magnitude/`shape`
0.002–0.024 at every e, against the 0.0695 bound. The joint construction stays size-pure at n = 80, as it was
at n = 1200.

**Orientation → `shape` and orientation → `delta` exceed the target statistic.** At e = 1.00 they are 0.464
and 0.354, against `angle` 0.300; at e = 0.50, 0.196 and 0.160 against 0.058. Under the default materiality
rule, the localization table labels both `projection_associated`: they first exceed the rule in the PLS latent
space (`shape` from e = 0.75, normalized excess 0.055–0.077; `delta` at e = 1.00, 0.062). That label must not be
carried into the paper. The
[magnitude-construction diagnostic](magnitude-construction-diagnostic-2026-09-10.md) (§1, §4) refuted the
projection-associated framing for orientation's cross-talk at n = 1200: the response is present within each omic
block, and under the recalibrated (null-dispersion) rule no pair localizes as projection-associated. That rule
is opt-in and was not used here. The orientation construction is a per-omic feature relocation, and its
`shape`/`delta` responses are reported as construction cross-talk. They are descriptive gate roles, not
acceptance targets.

**Shape's off-diagonals are construction-present.** shape → `angle` (0.138 at e = 1.00) is material in the
population-standardized geometry from e = 0.25 (normalized 0.236–0.436). shape → `delta` (0.330) is material
from e = 0.75 (0.054–0.065). Bending one vertex of a three-vertex path changes its direction and length, not
only its shape. Shape's construction moves `delta` more often than `shape`.

## 7. Drivers

Orientation attribution ran on the four nonzero orientation `power_primary` cells, observed component only
(`report/driver_report.csv`; contract echo in `report/report_contract.json`). For each transition (0→1, 1→2), at
e = 0.25 / 0.50 / 0.75 / 1.00:

| e | precision (0→1 / 1→2) | recall (0→1 / 1→2) | mean selected | within-replicate bootstrap sign stability |
|---|---|---|---|---|
| 0.25 | 0.996 / 0.996 | 0.375 / 0.379 | 20 | 0.787 / 0.786 |
| 0.50 | 1.000 / 1.000 | 0.194 / 0.193 | 20 | 0.813 / 0.810 |
| 0.75 | 1.000 / 1.000 | 0.131 / 0.132 | 20 | 0.835 / 0.835 |
| 1.00 | 1.000 / 1.000 | 0.101 / 0.102 | 20 | 0.855 / 0.856 |

The top-20 observed drivers are almost always true drivers. Recall falls with e because the truth set grows
(every relocated site and its cascade) while k stays at 20. `bootstrap_top_k_frequency_mean` is 0.040 at every
cell. That is top_k / n_features (20 / 498) by construction for a mean over all features, so it carries no
information here. The attribution figure is `report/phase4_attribution_stability.png` ("Within-replicate
bootstrap stability (observed component)").

> **No cross-replicate driver-stability claim is made.** Cross-replicate top-k Jaccard and sign agreement
> (`top_k_jaccard`, `sign_agreement` in `report/phase4_attribution.csv`) are `descriptive` only: every
> replicate index draws a fresh differential-indicator set, so the true driver set genuinely differs across
> replicates and matched seeds pair cells *within* a replicate index, not across them. A stability claim needs
> a design that holds the driver set fixed, which this study is not.

Driver precision is high even where `angle` power is low. Attribution describes the observed directional
contrast whether or not the `angle` test rejects, so in the case study a driver list is descriptive and is not
evidence of an orientation difference.

## 8. Construction limitations and case-study interpretability

- **This is a structural analogue, not an emulator.** InterSIM blocks are Gaussian/M-value with 367 + 131
  measured features against 80 samples. SEA-AD RNA and ATAC pseudobulks after VST are only approximately
  Gaussian and will have thousands of filtered features. The operating characteristics transfer as statements
  about the design (n, cell sizes, three stages, two blocks, stage-supervised PLS), not about RNA/ATAC
  distributions. The case-study repo must report its feature filtering beside any result.
- **Orientation power is eigengap- and rank-conditional.** At ρ = 0 the three stage means form a
  near-isotropic triangle, and its eigengap is small (anchor median 0.111). This is the generator's
  independent-baseline assumption, copied from Phase 5, not a consequence of three stages. The operating
  characteristics here are not measured at other baseline geometries, and Phase 5 found the orientation
  construction's realized contrast is not ρ-invariant. In the low eigengap tercile, `angle` power is at
  most 0.114 even at the construction maximum. At the construction maximum, power depends on whether CV retains more than three
  components, which a six-donor change in the cohort moves from 51 % to 8 % of replicates. The Phase 5 exit
  condition stands: the cohort's own recorded eigengap, not n, is what carries an orientation statement.
- **The orientation and shape constructions are bounded.** e = 1.00 relocates every differential site, so no
  larger construction exists. The HOLD is a statement that these constructions are not detectable at n = 80,
  not that the effect axis was too short.
- **The statistics are not mode-specific at n = 80.** The orientation construction is detected more often by
  `shape` and `delta` than by `angle`, and the shape construction more often by `delta` than by `shape`. Only
  magnitude is specific: the joint size change leaves `angle` and `shape` at the floor.
- **Deviation from the predeclared targets.** The orientation/`angle` and shape/`shape` power floors (0.80) are
  missed. Per design D6 this revises the **claim**, not the method or the Monte Carlo size. MOTCO's orientation
  and shape tests are not claimed to have useful power at this design.

**Per-statistic interpretability for the SEA-AD cohort (n = 80, ≥30 nuclei, RNA + ATAC):**

| statistic | rejection may be reported as evidence? | non-rejection may be reported as evidence? |
|---|---|---|
| `delta` | **Yes**, that the two sexes' ADNC trajectories differ. Type I is controlled (0.016–0.034). A rejection is consistent with a size difference, but it does not by itself exclude a direction or shape difference, since those constructions also reject `delta` at up to 0.35. | **Yes, against large size differences only.** Power is 0.956 for a 25 % joint size change and 0.542 for 10 %. A non-rejection argues against a size difference of about 25 % or more, and is uninformative below that. |
| `angle` | **Yes**, that the trajectories differ (Type I 0.024–0.040), but not specifically in direction. Orientation constructions reject `shape` and `delta` more often, and shape constructions reject `angle` at up to 0.14. Report it with the cohort's recorded eigengap. If that eigengap falls in the simulated low tercile (below about 0.13), a rejection is still valid at level α, but the test had little chance of one. | **No.** Power is at most 0.300, at the strongest construction the generator can build, and 0.03–0.09 below it. A non-rejection is uninformative at n = 80. |
| `shape` | **Yes**, that the trajectories differ (Type I 0.018–0.020), but not specifically in shape. The orientation construction rejects `shape` at 0.464, more than either `angle` or the shape construction's own 0.228. | **No.** Power is at most 0.228 at the construction maximum. A non-rejection is uninformative at n = 80. |

Using the ≥50-nuclei cohort (n = 74) leaves these statements unchanged for `delta` and makes `angle` and `shape`
weaker still. Adding a third standardized block would raise power for all three statistics, but the case study
has two blocks.

## 9. Reproduction

Reproduction never passes a worker-count override: `n_jobs` is part of every cell's parameter signature and
both configs' contracts forbid it (the runner exits 2). Parallelize across shards. Avoid nodes under memory
pressure: a stalled shard can be cancelled and resubmitted, and the signature-guarded resume skips completed
units.

```bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

# SLURM (partition/resource flags are cluster-specific); STUDY_N_JOBS must stay unset.
for P in magnitude:phase6-small-n-magnitude-2026-10-09 study:phase6-small-n-2026-10-09; do
  CFG=examples/trajectory_power_study/phase6_small_n_${P%%:*}.json; RUN=results/${P#*:}
  sbatch -p <partition> --cpus-per-task=1 --mem=2G --time=8:00:00 --array=0-99 \
      --export=ALL,STUDY_CONFIG=$(pwd)/$CFG,STUDY_OUT=$(pwd)/$RUN,N_SHARDS=100 \
      scripts/motco_study_array.sbatch
done

# After completion, per profile:
uv run python scripts/motco_study.py merge  --out-dir $RUN
uv run python scripts/motco_study.py report --config $CFG --out-dir $RUN
```

Each run commits `report/` and `PROVENANCE.txt` and leaves the shard and merged JSONL gitignored. The pilot that
fixed the effect axes is `results/phase6-small-n-pilot-2026-10-09/` (`NOTES.md`, `PROVENANCE.txt`).
