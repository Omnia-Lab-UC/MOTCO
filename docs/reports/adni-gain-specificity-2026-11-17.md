# ADNI gain and specificity — report (2026-11-17)

**Change:** `adni-gain-and-specificity` (MOTCO 0.25.0, no package change). **Reads:**
`results/adni-gain-specificity-2026-11-17/` (`adni_gain_specificity_decision.json`; per profile `settings.json`,
`preflight.json`, `rate_table.csv`, `paired_table.csv`, `cell_table.csv`, `accounting.json`;
`plants/reproduction.json`; `PROVENANCE.txt`). **Parents:** the [ADNI real-noise anchor](adni-real-noise-anchor-2026-11-17.md)
(its rate table and decision, bound by sha256) and the [operating characteristics at the ADNI design](adni-design-2026-11-02.md)
(its settings and per-unit p-values). **Concept page:**
[measurement visibility](../measurement-visibility.md#gain-and-specificity-at-adnis-design).

The ADNI anchor (item 21) found that `predictive_orthonormal` holds size on ADNI's blocks. Its decision fixed what else
reopening the option at a design needs: a gain over production there in simulation, and specificity there, each under
a decision declared before its reading. This change reads both, adds protected power, and reads what production itself
detects on ADNI's noise. It plants stretches, rotations and bends of the stage program on ADNI's own blocks, under
within-stage-permuted labels, and in the simulation of ADNI's design. One decision, fixed before any unit ran, reads
both.

**In one paragraph.** **The verdict is `power_loss`.**

- **The candidate loses target power on one plant.** On a 90° rotation of the analysed sample's own stage
  configuration, production's `angle` rejects in 0.285 of units at the full cohort and `predictive_orthonormal` in
  0.05. Production alone rejects 47 units and the candidate alone none, below S2's ceiling. At 750 participants the
  rates are 0.235 and 0.185.
- **It also fails specificity in simulation.** A doubled stretch is rejected by `shape` in 0.295 of units against
  production's 0.06, and two more off-target readings are below 0.025. On ADNI's in-sample plants it adds no off-target
  rejection, and it removes most of production's own rotation-into-`shape` leak.
- **Size holds and the gain is real.** Both arms hold size on every added control (S1). In simulation the candidate
  gains at the Bonferroni threshold in 5 of 8 target cells (G): `delta` 0.12 → 0.485 at a stretch of 0.5, `angle`
  0.075 → 0.325 at a rotation of 0.5.
- **The consequence fixed with the rule.** The reopening condition is not met for ADNI. `predictive_orthonormal` stays
  an opt-in diagnostic, and open item 1 closes for ADNI.
- **Production's statement** (on its own rates, whatever the verdict):
  - `delta` reaches the 0.80 floor from a stretch of e = 1 at 750 participants (0.995 in sample, 0.855 held out) and
    from e = 0.5 in sample at the full cohort (0.825).
  - `angle` and `shape` stay below the floor at every effect read; their largest rate is 0.40.
  - Its off-target family `leaks`: 10 of 32 readings are above the bound and 7 above the ceiling, mostly in-sample
    rotations read as `shape`.

## 1. Configuration and provenance

- **Export.** `adni_sexdiff/results/motco_input`, read in place through `ADNI_CONTRACT` and never copied between
  machines: 1,501 participants, the `p180` (133) and `nmr` (165) metabolite blocks. The workstation and `ing` each hold
  a byte-identical copy. The five contract files' hashes are in `plants/settings.json` and were unchanged after the
  run.
- **Profile `plants`: ADNI's blocks with a planted difference.** Each unit permutes the sex labels within stage, leaves
  blocks and covariates on their rows, and adds a known shift to the second group's rows, one vector per stage. No
  analysis ran under the real labels.

  | point | participants | Female CN / MCI / AD | Male CN / MCI / AD | held out |
  |---|---|---|---|---|
  | `n750` | 750 | 118 / 155 / 63 | 111 / 224 / 79 | the other 751 |
  | `cohort` | 1,501 | 237 / 310 / 125 | 223 / 448 / 158 | none |

  Units are the anchor driver's (its subsample, label and RRPP seeds, base seed 2039). The cohort anchor's 500 units
  are therefore item 21's, and every cell of a point at one replicate shares its sample, labels and permutation seed.
- **The plant.** A configuration is the three stage coefficients of `features ~ stage indicators + covariates`, in the
  analysis sample's pooled-standardized coordinates, with the covariates coded as `analyze_trajectories` codes them.
  `edit_stage_configuration` edits it:
  - a stretch about the centroid by 1 + e, whose intended difference is `delta`;
  - a rigid rotation by 90°·e, intended `angle`;
  - a bend of the middle stage by e mean steps at fixed path length, intended `shape`;
  - a translation by e mean steps along the net displacement, a group main effect with no interaction.

  The shift is the edited configuration minus the original, mapped back through the pooled SD.
- **Two references for a plant.** `sample` edits the analysis sample's own configuration. Under the label permutation
  both groups' expected configuration is the sample's, so the off-target truth is exact; the fitted span holds the
  plant by construction. `heldout` edits the configuration of the 751 participants outside the sample, a direction no
  fit was built on, available at `n750` only. Its off-target truth is not exact (§9.2).
- **Cells** (`<reference>_<mode>_e<effect>`; 500 replicates on a control, 200 on a plant, 999 permutations; 4,200
  units at `n750` and 3,000 at the cohort):

  | cells | reference | effects | points |
  |---|---|---|---|
  | `anchor`, `translation` (1.0) | —, `sample` | — | both |
  | stretch, rotation | `sample` | 0.25, 0.5, 1.0 | both |
  | bend | `sample` | 0.5, 1.0 | both |
  | stretch | `heldout` | 0.5, 1.0, 2.0 | `n750` |
  | rotation | `heldout` | 0.5, 1.0 | `n750` |
  | bend | `heldout` | 1.0 | `n750` |
  | out-of-span rotation (recorded, not judged) | `sample` | 0.25, 0.5 | both |

- **Profile `simulation`: stage-span plants at ADNI's design.** It uses the ADNI-design run's construction unchanged,
  at its calibrated strengths and ADNI's table. The generator's stage-span construction under the planted scale
  replaces the feature surgeries, so each plant's truth records its intended `delta`, `angle` and `shape`. Base seed
  2040. 2,600 units: the anchor and a translation (1.0) at 500, stretches at 0.5, 1.0 and 2.0, rotations at 0.25, 0.5
  and 1.0, and bends at 1.0 and 2.0, each at 200.
- **Arms.** Each arm is one `analyze_trajectories` call at the case study's configuration (five covariates in six
  design columns, rank floor 2, net displacement, the harness CV defaults) with the `df` residual correction, all at
  the unit's RRPP seed.
  - `pls_cv` cross-validates on the unit's data, planted or not; the cohort anchor is forced at rank 8.
  - `predictive_orthonormal` is forced at the same rank, so it shares the fit.
  - Both again under within-stage draws (`__within_stage`), recorded only.
- **Code and host.** Revision `2ecfe16`, pushed before any unit ran, from clean detached worktrees on both machines.
  Python 3.11.15 on the workstation and 3.11.16 on `ing`, with numpy 2.3.5, scipy 1.16.3, scikit-learn 1.8.0 and
  pandas 2.3.3 on both, BLAS pinned to one thread.
- **Where it ran.** The simulation ran on `ing` as array 944797: 100 tasks, all `COMPLETED`, 61.13 recorded
  core-hours against a cap of 90. The plants started on the workstation: 1,157 units, 7.46 core-hours. They then moved
  to `ing`, which holds the same export, and the job resumed the other 6,043 units under the same settings hash. In
  total 40.27 core-hours against a raised cap of 85. Only this run's own records were copied between machines
  (`PROVENANCE.txt`). The design had put the plants on the workstation alone, on the premise that the export existed
  only there.

## 2. Accounting, checks, pre-flights and the reproduction

- **Units.** Plants 7,200 of 7,200 and simulation 2,600 of 2,600 complete: no failed unit, no failed arm, no
  inconsistent unit, no duplicate. On ADNI's blocks no subsample and no label draw needed a redraw.
- **Checks on every unit.** Each forced arm reported its rank, every arm at one rank had the same model, each
  predictive arm measured two dimensions, and each within-stage arm observed its counterpart's statistics. The realized
  cell table was the point's. On ADNI's blocks only the second group's stage means moved, by the declared shift
  (`plant`). In simulation the truth recorded the stage-span construction at the cell's mode and effect
  (`construction`).
- **Plants pre-flight: passed.** On cohort replicates 0 and 1 the cross-validated fit selected rank 8 at a held-out
  stage AUROC of 0.599514, and the call forced at rank 8 gave identical statistics and p-values.
- **Simulation pre-flight: passed.** Anchor replicates 0 and 1 of the ADNI-design run, through its own evaluator under
  its committed settings, gave its committed p-values exactly: 0.29 / 0.153 / 0.195 and 0.187 / 0.746 / 0.387.
- **Anchor reproduction: passed.** The cohort anchor's counts for all four arms equal item 21's committed table, 12 of
  12 readings (`pls_cv` 25 / 32 / 28, `predictive_orthonormal` 29 / 20 / 18, and the two within-stage arms 35 / 38 /
  28 and 27 / 28 / 21 of 500).

## 3. The decision

Rule version 1 (`design.md` D12–D14). Rejection is p ≤ 0.05.

- **S1 (size).** Judged for both arms on the added controls: the `n750` anchor and the translation at both points (9
  readings), and the simulated anchor and translation (6). A family holds with at most one count above 34 of 500 and
  none above the ceiling (38 and 37).
- **S2, S3 and G.** One-sided exact paired tests of `predictive_orthonormal` against `pls_cv` on the units' discordant
  pairs.
  - S2 (protected power): the target statistic of all 30 plants. It holds with at most two loss p-values below 0.025
    and none below 0.05 / 30.
  - S3 (specificity): the off-target statistics of the 8 simulated plants (16 readings) and of the 16 in-sample plants
    (32). Each family holds with at most two excess p-values below 0.025 and none below 0.05 / m.
  - G (gain): at least one of the 8 simulated target readings gains at p ≤ 0.05 / 8.
- **Verdict.** `gate_failed`, then `production_excess` (production fails S1), `size_excess`, `power_loss` (S2),
  `specificity_excess` (S3), `no_gain` (G), else `reopened`. All failing criteria are listed.

| criterion | family | readings | observed | holds |
|---|---|---|---|---|
| S1, `pls_cv` | plants; simulation | 9; 6 | no count above the bound in either | yes |
| S1, `predictive_orthonormal` | plants; simulation | 9; 6 | no count above the bound in either | yes |
| S2 | all plants | 30 | 2 p-values below 0.025 (2 allowed), 1 below the ceiling 0.00167 | **no** |
| S3 | simulated plants | 16 | 3 below 0.025 (2 allowed), 1 below the ceiling 0.003125 | **no** |
| S3 | in-sample plants | 32 | none below 0.025 | yes |
| G | simulated targets | 8 | 5 gains at p ≤ 0.00625 (one needed) | yes |

**The verdict is `power_loss`**, with S2 and S3 failing. The gate passed and the decision is at paper grade (not
`indicative`). The consequence fixed with the rule: "predictive_orthonormal loses target power against production at
the named plants. The condition is not met for ADNI, and it stays an opt-in diagnostic. The measurement-visibility page
records which differences its two measured dimensions give up at ADNI's design. In every verdict the SEA-AD anchor's
reading stands for SEA-AD's design, keep_production stands, and the out-of-span cells are reported beside the verdict
and were not judged."

**What the verdict licenses** (written into the decision beside item 21's scoping rule S1–S5):

- **L1.** "The reading speaks for ADNI's export, the two design points and the simulated table. Nothing is extended
  to 80, 200 or 500 participants, where only size was read."
- **L2.** "At the full cohort, power and specificity were read on in-sample plants and in simulation. No held-out
  plant exists there. The held-out readings are those of 750 participants."
- **L3.** "reopened permits a proposal. It changes no default and adds no option. The later change carries the rank
  rule open item 1 asks for." (Not reached under this verdict.)
- **L4.** "Cohorts do not vouch for one another. The SEA-AD veto stands for SEA-AD's design, and ROSMAP has no
  reading."
- **L5.** "The out-of-span cells license nothing. Outside the stage span the split test's complement is the detector,
  and its reading on ADNI's noise is not this change's."

## 4. Size on the added controls (S1)

Rejections of 500 for `delta` / `angle` / `shape`:

| control | `pls_cv` | `predictive_orthonormal` |
|---|---|---|
| `n750` anchor | 12 / 13 / 15 | 15 / 13 / 14 |
| `n750` translation | 11 / 17 / 14 | 10 / 22 / 13 |
| cohort translation | 24 / 30 / 23 | 27 / 18 / 25 |
| simulated anchor | 21 / 33 / 22 | 24 / 19 / 28 |
| simulated translation | 21 / 33 / 21 | 22 / 16 / 26 |

No count is above the bound of 34. **A group main effect along the stage program does not reach the interaction
statistics at ADNI's design**, in either measurement, on real noise or in simulation. The `n750` anchor, a point item
21 did not read, has every rate between 0.024 and 0.030.

## 5. Protected power (S2)

Target rates of 200 units, production then the candidate, with the discordant units: b = production alone rejects, c =
the candidate alone. p is the one-sided exact loss p-value.

| profile, point | cell | statistic | `pls_cv` | `predictive_orthonormal` | b / c | p (loss) |
|---|---|---|---|---|---|---|
| simulation | stretch 0.5 / 1.0 / 2.0 | `delta` | 0.12 / 0.535 / 0.96 | 0.485 / 0.99 / 1.0 | 5 / 78, 0 / 91, 0 / 8 | 1.0 |
| simulation | rotation 0.25 / 0.5 / 1.0 | `angle` | 0.07 / 0.075 / 0.065 | 0.13 / 0.325 / 0.40 | 13 / 25, 7 / 57, 4 / 71 | 0.98, 1.0, 1.0 |
| simulation | bend 1.0 / 2.0 | `shape` | 0.045 / 0.04 | 0.08 / 0.105 | 8 / 15, 8 / 21 | 0.95, 0.996 |
| `n750`, held out | stretch 0.5 / 1.0 / 2.0 | `delta` | 0.22 / 0.855 / 1.0 | 0.24 / 0.935 / 1.0 | 16 / 20, 8 / 24, 0 / 0 | 0.80, 0.999, 1.0 |
| `n750`, held out | rotation 0.5 / 1.0 | `angle` | 0.075 / 0.345 | 0.20 / 0.465 | 8 / 33, 20 / 44 | 1.0, 0.999 |
| `n750`, held out | bend 1.0 | `shape` | 0.055 | 0.055 | 6 / 6 | 0.61 |
| `n750`, in sample | stretch 0.25 / 0.5 / 1.0 | `delta` | 0.175 / 0.665 / 0.995 | 0.275 / 0.89 / 1.0 | 13 / 33, 3 / 48, 0 / 1 | 0.999, 1.0, 1.0 |
| `n750`, in sample | rotation 0.25 / 0.5 / **1.0** | `angle` | 0.025 / 0.14 / **0.235** | 0.23 / 0.305 / **0.185** | 2 / 43, 3 / 36, **14 / 4** | 1.0, 1.0, **0.015** |
| `n750`, in sample | bend 0.5 / 1.0 | `shape` | 0.115 / 0.22 | 0.14 / 0.27 | 10 / 15, 7 / 17 | 0.89, 0.99 |
| cohort, in sample | stretch 0.25 / 0.5 / 1.0 | `delta` | 0.255 / 0.825 / 1.0 | 0.355 / 0.915 / 1.0 | 18 / 38, 2 / 20, 0 / 0 | 0.998, 1.0, 1.0 |
| cohort, in sample | rotation 0.25 / 0.5 / **1.0** | `angle` | 0.10 / 0.23 / **0.285** | 0.36 / 0.345 / **0.05** | 8 / 60, 6 / 29, **47 / 0** | 1.0, 1.0, **7.1e-15** |
| cohort, in sample | bend 0.5 / 1.0 | `shape` | 0.185 / 0.40 | 0.245 / 0.565 | 14 / 26, 11 / 44 | 0.98, 1.0 |

- **Two readings are below 0.025, both the in-sample rotation of 90°.** At the cohort the loss is complete (47 units
  against none) and below the ceiling of 0.05 / 30, so S2 fails. At `n750` it is 14 against 4.
- **Everywhere else the candidate reads at least as often as production.** That includes the in-sample rotations of
  22.5° and 45° (0.025 → 0.23 at `n750`, 0.10 → 0.36 at the cohort) and the held-out rotations (0.075 → 0.20 at 45°).
- **What distinguishes the failing plant.** A 90° rotation of one group cancels much of the pooled stage signal. Under
  it the CV selects rank 2 in 143 of 200 units at `n750` and 115 at the cohort, at a median held-out AUROC of 0.515
  and 0.518 (§9.6). Among the cohort's 85 units with a rank above 2, production rejects `angle` in 51 and the candidate
  in 6; among the 115 at rank 2, in 6 and 4 (§9.5). The decision explains no failure. This split is descriptive.

## 6. Specificity (S3)

**Simulated plants: fails.** Off-target rates of 200 units, production then the candidate, with the one-sided exact
excess p-value:

| cell | statistic | `pls_cv` | `predictive_orthonormal` | b / c | p (excess) |
|---|---|---|---|---|---|
| stretch 0.5 | `angle` / `shape` | 0.045 / 0.045 | 0.035 / 0.06 | 9 / 7, 8 / 11 | 0.77, 0.32 |
| stretch 1.0 | `angle` / `shape` | 0.04 / 0.045 | 0.055 / **0.11** | 8 / 11, **6 / 19** | 0.32, **0.0073** |
| stretch 2.0 | `angle` / `shape` | 0.025 / 0.06 | 0.035 / **0.295** | 5 / 7, **7 / 54** | 0.39, **2.2e-10** |
| rotation 0.25 | `delta` / `shape` | 0.03 / 0.045 | 0.05 / 0.04 | 3 / 7, 9 / 8 | 0.17, 0.69 |
| rotation 0.5 | `delta` / `shape` | 0.025 / 0.06 | 0.06 / 0.05 | 3 / 10, 12 / 10 | 0.046, 0.74 |
| rotation 1.0 | `delta` / `shape` | 0.035 / 0.07 | **0.08** / 0.095 | **2 / 11**, 11 / 16 | **0.011**, 0.22 |
| bend 1.0 | `delta` / `angle` | 0.03 / 0.08 | 0.06 / 0.05 | 3 / 9, 16 / 10 | 0.073, 0.92 |
| bend 2.0 | `delta` / `angle` | 0.03 / 0.085 | 0.055 / 0.055 | 3 / 8, 17 / 11 | 0.11, 0.91 |

Three readings are below 0.025, one more than allowed, and the doubled stretch's `shape` is below the ceiling of
0.05 / 16. **The candidate reads a stretch of the simulated program as a change of shape**, more so as the stretch
grows (0.06 → 0.11 → 0.295), where production stays at 0.045–0.06.

**In-sample plants on ADNI's blocks: holds.** No excess p-value is below 0.025 (the smallest is 0.125). Off-target
rejections of 200, production / candidate:

| cell | `n750` | cohort |
|---|---|---|
| stretch 0.25: `angle`, `shape` | 4 / 5, 7 / 12 | 16 / 12, 13 / 10 |
| stretch 0.5: `angle`, `shape` | 7 / 6, 12 / 12 | 19 / 15, 16 / 11 |
| stretch 1.0: `angle`, `shape` | 21 / 11, 16 / 21 | 31 / 18, 19 / 18 |
| rotation 0.25: `delta`, `shape` | 7 / 4, 10 / 12 | 8 / 10, 24 / 17 |
| rotation 0.5: `delta`, `shape` | 10 / 5, 41 / 15 | 6 / 9, 44 / 18 |
| rotation 1.0: `delta`, `shape` | 13 / 1, 79 / 20 | 3 / 0, 61 / 4 |
| bend 0.5: `delta`, `angle` | 7 / 7, 4 / 4 | 12 / 12, 12 / 12 |
| bend 1.0: `delta`, `angle` | 10 / 7, 5 / 6 | 10 / 10, 19 / 14 |

On ADNI's noise the candidate removes most of production's rotation-into-`shape` leak: at 90° production alone
rejects `shape` in 59 units at `n750` and 57 at the cohort, the candidate alone in none.

## 7. Gain (G) and the recorded gains

**G holds.** Five of the eight simulated target readings gain at p ≤ 0.00625 (§5): the three stretches (`delta`
0.12 → 0.485, 0.535 → 0.99, 0.96 → 1.0) and the rotations of 0.5 and 1.0 (`angle` 0.075 → 0.325 and 0.065 → 0.40).
The rotation of 0.25 and both bends do not. Pooled over the eight, the candidate alone rejects 366 units and
production alone 45, a mean rate difference of +0.20.

**Recorded beside G** (same rule, no criterion):

- **Held-out plants (6 readings): holds.** Three gain, the stretch of 1.0 and both rotations (p 0.0035, 5.6e-5,
  0.0018). Pooled, c = 127 against b = 58, +0.0575.
- **In-sample plants (16 readings): holds.** Eight gain, all on stretches, rotations up to 45° and the cohort's bend
  of 1.0. Pooled, c = 414 against b = 158, +0.08.

The gain the simulation shows is also on ADNI's noise. The verdict turns on what is lost, not on what is gained.

## 8. Production's statement

Read on `pls_cv` alone. It is written whatever the verdict and never changes it. **Power** (target rejection rates;
the floor is 0.80). A held-out plant bounds the detectable effect from above, an in-sample one from below:

| point, reference | `delta` (stretch) | `angle` (rotation) | `shape` (bend) |
|---|---|---|---|
| `n750`, in sample | 0.175 / 0.665 / **0.995** at 0.25 / 0.5 / 1.0: from e = 1 | 0.025 / 0.14 / 0.235: none | 0.115 / 0.22 at 0.5 / 1.0: none |
| `n750`, held out | 0.22 / **0.855** / 1.0 at 0.5 / 1.0 / 2.0: from e = 1 | 0.075 / 0.345 at 0.5 / 1.0: none | 0.055 at 1.0: none |
| cohort, in sample | 0.255 / **0.825** / 1.0 at 0.25 / 0.5 / 1.0: from e = 0.5 | 0.10 / 0.23 / 0.285: none | 0.185 / 0.40 at 0.5 / 1.0: none |

- **`delta`.** At 750 participants both references reach the floor at a stretch of e = 1, so the detectable stretch
  of ADNI's adjusted stage program is bracketed at e = 1 on this ladder. At the cohort the in-sample plant reaches it
  at e = 0.5.
- **`angle` and `shape`.** Neither reaches the floor on either reference at either point up to e = 1. On ADNI's
  noise a non-rejection of `angle` or `shape` stays uninformative, as item 9 found in simulation.

**Off-target: `leaks`.** Production's 32 off-target readings on the in-sample plants, as one family under the
count-plus-ceiling rule at 200 units (at most two above 16 of 200, none above 20), have 10 above the bound and 7
above the ceiling:

- `shape` under in-sample rotations: 41 and 79 of 200 at `n750` (e = 0.5, 1.0); 24, 44 and 61 at the cohort.
- `angle` under stretches: 21 at `n750` (e = 1.0); 19 and 31 at the cohort (e = 0.5, 1.0).
- The cohort's stretch of 1.0 by `shape` (19) and bend of 1.0 by `angle` (19).

**On ADNI's noise, a production rejection names that the trajectories differ, not how.** A rotation of the stage
program reaches `shape` more often than `angle` (0.395 against 0.235 at `n750`, e = 1).

## 9. Descriptive readings

None changes the verdict or production's statement.

### 9.1 The bracket

Target rejections of 200 at `n750`, in sample against held out, at the same mode and effect:

| plant | `pls_cv` | `predictive_orthonormal` |
|---|---|---|
| stretch 0.5 (`delta`) | 133 / 44 | 178 / 48 |
| stretch 1.0 (`delta`) | 199 / 171 | 200 / 187 |
| rotation 0.5 (`angle`) | 28 / 15 | 61 / 40 |
| rotation 1.0 (`angle`) | 47 / 69 | 37 / 93 |
| bend 1.0 (`shape`) | 44 / 11 | 54 / 11 |

The in-sample plant is read more often, except the rotation of 90°: there the held-out plant is read more often by
both arms. Rotating a configuration that is mostly the reference's noise does not cancel the sample's own pooled
stage signal.

### 9.2 Held-out off-target rates and the reference's reliability

The held-out configuration's split-half correlation over the 500 `n750` anchor units has a median of 0.166, and its
Spearman–Brown step-up to the full reference a median of 0.285. About 0.3 of a held-out plant's energy is the
cohort's stage program and the rest is estimation noise, as design D5 stated. The second eigenvalue of the halves'
cross-product has a median ratio to the first of 0.000. Off-target rejections of 200, production / candidate:

| held-out cell | off-target |
|---|---|
| stretch 0.5 | `angle` 8 / 9, `shape` 9 / 5 |
| stretch 1.0 | `angle` 24 / 14, `shape` 22 / 18 |
| stretch 2.0 | `angle` 81 / 38, `shape` 44 / 53 |
| rotation 0.5 | `delta` 39 / 69, `shape` 24 / 35 |
| rotation 1.0 | `delta` 163 / 176, `shape` 64 / 91 |
| bend 1.0 | `delta` 7 / 4, `angle` 4 / 3 |

A held-out rotation is read by `delta` in up to 176 of 200 units. This is the plant's own noise: an edit of a noisy
configuration carries an unintended size change. None of these rates was judged.

### 9.3 The out-of-span cells (recorded, not judged)

A rotation of the net displacement by 90°·e toward a seeded direction orthogonal to the sample's stage span.
Rejections of 200 for `delta` / `angle` / `shape`, and the discordant units b / c per statistic:

| point, effect | `pls_cv` | `predictive_orthonormal` | b / c (`delta`; `angle`; `shape`) | median CV rank | plant in span / in predictive subspace |
|---|---|---|---|---|---|
| `n750`, 0.25 | 28 / 24 / 22 | 66 / 77 / 147 | 8 / 46; 15 / 68; 1 / 126 | 19 | 0.466 / 0.374 |
| `n750`, 0.5 | 152 / 167 / 108 | 180 / 139 / 196 | 6 / 34; 48 / 20; 0 / 88 | 19 | 0.611 / 0.436 |
| cohort, 0.25 | 56 / 127 / 43 | 64 / 138 / 193 | 34 / 42; 37 / 48; 0 / 150 | 19 | 0.523 / 0.430 |
| cohort, 0.5 | 199 / 200 / 117 | 193 / 192 / 200 | 7 / 1; 8 / 0; 0 / 83 | 19 | 0.637 / 0.463 |

- **On ADNI's noise the fit follows the plant.** The CV selects rank 19 (the cap), and the fitted span keeps about
  half the plant.
- **Both arms reject; the candidate reads the rotation as `shape` above all.** Its `shape` rate is 0.735–1.0, against
  production's 0.11–0.585.
- **The split test's complement is the designated detector outside the stage span** (L5). Its reading on ADNI's noise
  is not this change's.

### 9.4 The within-stage arms

- **Size.** On the added controls `pls_cv__within_stage` has one count above the bound: the cohort translation's
  `angle`, 35 of 500, within its family's allowance. Its simulated family has none. `predictive_orthonormal__within_stage`
  has none in either family. At the cohort anchor both reproduce item 21's counts (§2): 35 / 38 / 28 for production,
  whose `angle` item 21 found above that point's ceiling.
- **Against the unrestricted draws on the same units.** 20 of the 204 paired comparisons have a two-sided p-value below
  0.05. In every one the within-stage arm rejects more, by 6 to 77 units, with at most one the other way.
  - The largest are the in-sample rotations' `angle`: at the cohort's 90° production goes 0.285 → 0.67 and the
    candidate 0.05 → 0.19; at `n750`, 0.235 → 0.385 and 0.185 → 0.325.
  - Off-target, at `n750`'s 90° rotation `shape` goes 0.395 → 0.435 (production) and 0.10 → 0.155 (candidate).
  - In simulation the stretch of 1.0's `delta` goes 0.535 → 0.605 for production.

### 9.5 Rates by rank class

Rejections by `pls_cv` / `predictive_orthonormal` for `delta`, `angle`, `shape`, among units at CV rank 2 (= K − 1,
where the predictive subspace is the fitted span) and above it:

| cell | units at rank 2 | at rank 2 | units above | above rank 2 |
|---|---|---|---|---|
| `n750` anchor | 53 | 1 / 0, 1 / 2, 1 / 1 | 447 | 11 / 15, 12 / 11, 14 / 13 |
| `n750` rotation 0.5 | 67 | 3 / 2, 1 / 2, 9 / 4 | 133 | 7 / 3, 27 / 59, 32 / 11 |
| `n750` rotation 1.0 | 143 | 3 / 0, 14 / 13, 43 / 13 | 57 | 10 / 1, 33 / 24, 36 / 7 |
| cohort rotation 0.5 | 125 | 3 / 4, 12 / 9, 9 / 5 | 75 | 3 / 5, 34 / 60, 35 / 13 |
| cohort rotation 1.0 | 115 | 1 / 0, 6 / 4, 12 / 3 | 85 | 2 / 0, 51 / 6, 49 / 1 |
| simulated rotation 1.0 | 34 | 2 / 2, 1 / 3, 4 / 2 | 166 | 5 / 14, 12 / 77, 10 / 17 |

- **The cohort's 90° loss sits above rank 2.** Production's `angle` reads 51 of those 85 units and the candidate 6.
- **At 45° and in simulation the candidate's gain is also above rank 2**, where the two coordinate systems differ.

### 9.6 The fit under a plant, the plant's size and the scale drift

Median CV rank and held-out stage AUROC per cell (the cohort anchor is forced at rank 8, AUROC 0.5995 from the
pre-flight):

| cell | `n750` rank, AUROC | cohort rank, AUROC |
|---|---|---|
| anchor | 5, 0.578 | 8 (forced) |
| translation 1.0 | 5, 0.576 | 7, 0.599 |
| stretch 0.25 / 0.5 / 1.0 (in sample) | 7, 0.626 / 9, 0.674 / 13, 0.750 | 8, 0.634 / 8, 0.671 / 18, 0.745 |
| rotation 0.25 / 0.5 / 1.0 (in sample) | 5, 0.566 / 3.5, 0.549 / 2, 0.515 | 7, 0.590 / 2, 0.561 / 2, 0.518 |
| bend 0.5 / 1.0 (in sample) | 5, 0.562 / 5, 0.554 | 7, 0.584 / 7, 0.570 |
| stretch 0.5 / 1.0 / 2.0 (held out) | 7, 0.617 / 7, 0.668 / 9, 0.763 | — |
| rotation 0.5 / 1.0 (held out) | 6, 0.591 / 8, 0.626 | — |
| bend 1.0 (held out) | 5, 0.563 | — |
| out-of-span 0.25 / 0.5 | 19, 0.632 / 19, 0.696 | 19, 0.682 / 19, 0.730 |

- **A real difference moves a stage-supervised fit.** A stretch raises the rank and the AUROC. An in-sample rotation
  of 45° or more lowers them, to rank 2 in 143 of 200 units at `n750` (90°) and in 125 and 115 at the cohort (45°,
  90°).
- **Plant size.** A unit in-sample stretch's shift has a median norm of 1.066 / 0.876 / 1.010 pooled SDs at CN / MCI /
  AD at `n750` (largest entry 0.209). A held-out stretch's is 1.050 / 0.875 / 1.031, and a 90° in-sample rotation's
  1.508 / 1.239 / 1.429. The translation moves every stage by 1.611 at `n750` and 1.334 at the cohort.
- **Scale drift.** The largest relative change of a pooled SD (median unit) is 0.13 % to 1.6 % in every cell at
  effects up to 1.0, and 3.4 % for the held-out stretch of 2.0 (largest unit 6.0 %).

### 9.7 Capture

Share of energy inside the fitted span / inside the predictive subspace (medians):

| configuration or plant | `n750` | cohort |
|---|---|---|
| the sample's own configuration (anchor) | 0.925 / 0.548 | 0.957 / 0.396 |
| the held-out configuration (anchor) | 0.536 / 0.109 | — |
| in-sample stretch 0.25 / 1.0 | 0.962 / 0.428, 0.986 / 0.214 | 0.964 / 0.401, 0.988 / 0.152 |
| in-sample rotation 0.25 / 1.0 | 0.904 / 0.574, 0.600 / 0.539 | 0.961 / 0.465, 0.623 / 0.623 |
| in-sample bend 0.5 / 1.0 | 0.874 / 0.583, 0.792 / 0.569 | 0.917 / 0.460, 0.904 / 0.470 |
| held-out stretch 0.5 / 2.0 | 0.670 / 0.174, 0.886 / 0.224 | — |
| held-out rotation 0.5 / 1.0 | 0.547 / 0.125, 0.666 / 0.118 | — |
| held-out bend 1.0 | 0.329 / 0.089 | — |

The predictive subspace keeps 0.17–0.22 of a held-out stretch and reads it at least as often as the scores (§9.1): as
in the probes, the share is not a measure of what that subspace detects.

### 9.8 In simulation: intended against observed geometry

The CV selects rank 19 in 494 of 500 anchor units, at a median held-out AUROC of 0.608. The population path length is
0.67 in planted-standardized units. Medians of the observed statistics (`delta` / `angle` / `shape`) beside each
plant's intended difference:

| cell | intended | `pls_cv` | `predictive_orthonormal` |
|---|---|---|---|
| anchor | 0 / 0° / 0 | 0.428 / 30.2° / 0.123 | 0.055 / 4.1° / 0.072 |
| stretch 1.0 | `delta` 0.674 | 1.347 / 26.8° / 0.111 | 0.308 / 2.9° / 0.080 |
| stretch 2.0 | `delta` 1.347 | 2.885 / 26.5° / 0.113 | 0.675 / 2.6° / 0.091 |
| rotation 1.0 | `angle` 90° | 0.450 / 30.2° / 0.157 | 0.090 / 14.1° / 0.114 |
| bend 2.0 | `shape` 0.490 | 0.437 / 31.3° / 0.122 | 0.052 / 3.9° / 0.093 |

The two arms measure in different coordinates, so their statistics are not on one scale. Production's `angle` stays
at its null median under a 90° rotation; the candidate's rises from 4.1° to 14.1°.

## 10. What follows

- **The reopening condition is not met for ADNI.** `predictive_orthonormal` stays an opt-in diagnostic, and open item
  1 closes for ADNI. Open item 2 continues with ROSMAP for production alone.
- **What its two measured dimensions give up at ADNI's design** (recorded on the measurement-visibility page):
  - a rotation of the sample's own stage program by 90°, where the pooled stage signal cancels: production reads its
    `angle` in units whose CV rank is above 2, and the candidate barely (§9.5);
  - in simulation, specificity on stretches, which it reads as `shape` more often as they grow.
- **The gains are real on both kinds of noise.** They are recorded, and they do not outweigh a loss under the rule.
- **Production's statement is the case study's real-noise power statement.** `adni_sexdiff` may cite it: on ADNI's
  noise `delta` detects a stretch of the adjusted stage program from e = 1 at 750 participants (from e = 0.5 in
  sample at the cohort); `angle` and `shape` reach the floor nowhere; and a rejection names a difference, not its
  kind.
- **The out-of-span cells are the starting point of the split test's reading on ADNI's noise.**
- **Standing.** The SEA-AD anchor's reading stands for SEA-AD's design, and `keep_production` stands.

## 11. Limitations

- **Subsamples of one cohort are not independent cohorts.** Two `n750` subsamples share about half their
  participants, and every rate is conditional on ADNI's 1,501. Independent real noise needs ROSMAP.
- **The two references bracket; neither is the truth.** An in-sample plant lies where the stage-supervised fit already
  looks, and a held-out one is about a third signal (§9.2). The cohort has no held-out plant (L2).
- **The simulation is a structural analogue.** Its CV selects rank 19 against ADNI's 8, its blocks are InterSIM's, and
  its stage program is weak in pooled-SD units. Every failing specificity reading is simulated; on ADNI's in-sample
  plants the candidate holds.
- **A rotation and a bend are exact to the order of the scale drift.** That is under 1.6 % of a pooled SD at effects
  up to 1.0 (median unit), and both arms read the same data.
- **The plant moves the fit.** A cell's rank and AUROC are not the anchor's (§9.6). The failing rotation is the cell
  where the pooled stage signal cancels and the CV falls to rank 2 in most units.
- **The family rules have false-alarm rates.** At exact size a correctly sized arm fails the plants size family with
  probability at most 0.07 and the simulated one 0.06. The verdict does not rest on a reading near a boundary: both
  failing families have a reading far below their ceiling (7.1e-15 against 0.00167; 2.2e-10 against 0.003125).
- **One omic type.** ADNI's blocks are metabolites at 0.2 to 0.4 features per participant.
