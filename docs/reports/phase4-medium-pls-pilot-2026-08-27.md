# Phase 4 medium PLS pilot

> **Partially superseded (2026-08-28) — orientation `angle` values only.** The orientation estimator used for
> this pilot anchored PC1's sign on the centered first-stage row, which flips the sign for bent trajectories.
> Every `angle` figure below therefore carries that defect. The corrected re-run is
> [`orientation-sign-anchor-2026-08-28.md`](orientation-sign-anchor-2026-08-28.md); corrected orientation
> `angle` power is 0.64/0.64/0.67/0.68 across effects 0.25 → 1.00, against 0.59/0.59/0.64/0.65 here.
>
> **The pilot's conclusions are not overturned.** The shortfall against the 0.80 floor survives correction,
> the Type I and translation controls still hold, the magnitude and shape conclusions are untouched, and the
> realized-geometry figures (81.1° population, 57.5° latent) substantially stand. What changes is the
> *explanation* offered below: the "why orientation is flat" section's reasoning is superseded, since the
> corrected estimator leaves the near-180° null mass unchanged at 21/700.
>
> The original numbers in this report are left exactly as run. Nothing below has been edited.


**Run date:** 2026-08-27
**Configuration:** `examples/trajectory_power_study/phase4_pilot_100x199.json`
**Code revision:** `073bb2f` on `feat/phase4-medium-pls-pilot` (parent `fabdc18`). The pilot executed from
this working tree before it was committed; the per-file fingerprints taken at launch are in
`results/phase4-2026-08-27/PROVENANCE.txt`, and the report outputs regenerate byte-identically from this
revision.
**Environment:** Python 3.11.15, numpy 2.3.5, pandas 2.3.3, scikit-learn 1.8.0, scipy 1.16.3, motco 0.6.0; no R at runtime
**Supersedes for gate purposes:** [`mvalue-pls-pilot-2026-07-30.md`](mvalue-pls-pilot-2026-07-30.md), which predates the
corrected shape estimator, realized-geometry diagnostics, and orientation attribution. That report is retained
unchanged as historical evidence and is **not** used for the Phase 5 decision.

## Phase 5 gate decision: **HOLD**

> Mandatory gate failed: `mandatory_power[orientation,angle]` — top-effect power `0.650` against the
> predeclared floor of `0.800`. The curve itself is monotone within tolerance; the level is short, and the
> cause is **not yet identified** (see "Why orientation is flat").

Every other mandatory gate passed. Nothing here is a Monte Carlo accident: at 100 replicates the standard
error at a rate of 0.65 is 0.048, so the observed power is more than three standard errors below the floor.
Raising the replicate count would sharpen the estimate without moving it.

The full decision, every observation, and its Monte Carlo standard error are in
`results/phase4-2026-08-27/report/phase4_gate_decision.json`.

## Run configuration

| Parameter | Value |
|---|---:|
| Integration | Pooled PLS, M-value methylation |
| PLS cross-validation | `cv1_splits=3`, `cv2_splits=4`, `n_repeats=5`, `max_components=20`, `random_state=1203` |
| Samples per replicate | 300 |
| Stages | 4 |
| Replicates per cell | 100 |
| Permutations per test | 199 |
| Effect sizes | 0.00, 0.25, 0.50, 0.75, 1.00 |
| Trajectory modes | magnitude, orientation, shape, translation |
| Seed policy | matched primary seeds, one shared zero-effect anchor |
| Attribution | all nonzero primary orientation cells; 100 bootstraps, `top_k=20` |
| Significance level | 0.05 |
| Cells | 19 |
| Total work units | 1,900 |
| Failed work units | 0 |
| Compute | 23.3 core-hours (attribution 0.45% of it) |

All 1,900 units completed, every parameter signature matched enumeration, no records were missing or
duplicated, and every completed record carried selected-component and realized-geometry metadata. All 400
eligible orientation replicates produced valid attribution diagnostics; none failed.

## Operating characteristics

Rejection rates at alpha 0.05, 100 replicates per entry. Bold is the statistic each mode is constructed to
move. Every mode's `0.00` row is the **same** shared zero-effect anchor cell, so those four rows are one
measurement, not four.

| Mode | Effect | Delta | Angle | Shape |
|---|---:|---:|---:|---:|
| magnitude | 0.25 | **1.00** | 0.04 | 0.04 |
| magnitude | 0.50 | **1.00** | 0.03 | 0.01 |
| magnitude | 0.75 | **1.00** | 0.01 | 0.00 |
| magnitude | 1.00 | **1.00** | 0.01 | 0.00 |
| orientation | 0.25 | 0.60 | **0.59** | 0.97 |
| orientation | 0.50 | 0.69 | **0.59** | 1.00 |
| orientation | 0.75 | 0.67 | **0.64** | 0.99 |
| orientation | 1.00 | 0.65 | **0.65** | 0.99 |
| shape | 0.25 | 0.80 | 0.51 | **0.94** |
| shape | 0.50 | 0.97 | 0.69 | **1.00** |
| shape | 0.75 | 1.00 | 0.77 | **1.00** |
| shape | 1.00 | 1.00 | 0.83 | **1.00** |
| translation | 0.25 | 0.03 | 0.03 | 0.06 |
| translation | 0.50 | 0.08 | 0.07 | 0.07 |
| translation | 0.75 | 0.07 | 0.05 | 0.04 |
| translation | 1.00 | 0.07 | 0.05 | 0.04 |
| *(shared anchor)* | 0.00 | 0.08 | 0.04 | 0.03 |

Source: `report/phase4_operating.csv`.

### Type I control

Eighteen one-sided control tests were predeclared: the `none` baseline plus every translation effect level,
each contributing one test per statistic. **All passed.** The one-sided bound at `n = 100` is
`0.05 + 2·sqrt(0.05·0.95/100) = 0.0936`; the largest observed control rate anywhere was `0.09`
(translation negative control, shape). No exceedance occurred, so the marginal-exceedance confirmation rule
was not engaged.

This is the strongest result in the pilot. Translation — a pure location offset — leaves all three
statistics at nominal level at every effect size, which is what a correctly calibrated test must do.

### Magnitude

`delta` reaches 1.00 at every nonzero effect and is monotone within tolerance. Both mandatory off-diagonal
controls pass: magnitude's `angle` (0.01) and `shape` (0.00) at effect 1.00 sit *below* alpha, not merely
below the inflation bound. The M-value contract established in July holds under the corrected estimator.

### Shape

`shape` power is 0.94–1.00 across the sweep and monotone. The corrected pairwise Procrustes estimator
detects the constructed interior bend easily, so the detectability requirement the gate encodes is met with
room to spare.

### Orientation — the failing gate

`angle` power is 0.59, 0.59, 0.64, 0.65 across effects 0.25 → 1.00. It is monotone within tolerance but
essentially **flat**: quadrupling the requested effect buys about six points of power. A construction whose
power does not respond to its own effect size is not simply underpowered — something other than the
requested signal is governing the test. The next section shows the obvious candidate (a sampling noise
floor) is contradicted by the records.

## Why orientation is flat: what the geometry does and does not explain

Realized geometry at effect 1.00, joint scope (`report/phase4_geometry.csv`). Angles are in degrees;
`delta` and `shape` are only comparable **within** a column, never across columns.

| Mode | Statistic | Population (standardized) | Observed (standardized) | PLS latent |
|---|---|---:|---:|---:|
| magnitude | angle | 39.9 | 45.8 | 32.0 |
| magnitude | delta | 27.5 | 27.2 | 26.8 |
| orientation | angle | 81.1 | 81.0 | 57.5 |
| orientation | delta | 1.69 | 1.69 | 2.19 |
| shape | angle | 85.0 | 85.8 | 65.2 |
| shape | delta | 8.96 | 8.89 | 9.95 |
| translation | angle | 0.0 | **16.8** | 8.9 |
| translation | delta | 0.0 | 0.38 | 0.37 |

The construction is strong and survives integration: 81° at the population level, still 57° in the latent
space, against a null floor of 8.9° (translation, whose population angle is exactly 0°). On a simple
detection-threshold reading that separation should yield near-complete power, and it does not.

**The obvious explanation is wrong, and the records say so.** Within the orientation cell at effect 1.00,
the replicates that *fail* to reject have a **larger** mean observed latent angle (66.7°) than the ones that
reject (52.5°). A noise-floor or underpowered-signal account predicts the opposite ordering. The same
inversion appears in shape at effect 1.00 (83.3° non-rejecting versus 61.5° rejecting), while the
translation null behaves normally (30.7° rejecting versus 7.8° non-rejecting) — so the inversion is specific
to cells that carry real signal, not a general property of the statistic.

Two candidate explanations are ruled out. Latent dimensionality is not responsible: splitting the cell by
selected components gives an angle rejection rate of 0.71 at `lv = 2` (n = 7) against 0.65 at `lv = 3`
(n = 93). And it is not Monte Carlo noise: the effect is a 14-degree gap across 65 versus 35 replicates.

What remains is that the RRPP permutation null for `angle` **co-varies with the observed statistic within a
replicate** — that is, the statistic is not pivotal across replicates, so a replicate whose latent geometry
inflates the observed angle inflates its own null by at least as much. The pilot cannot confirm this because
`include_null_distributions` was off, so the per-replicate null quantiles were not persisted. Confirming or
refuting it is the first item of pre-Phase-5 work; see
[Phase 5 readiness](https://github.com/Omnia-Lab-UC/MOTCO/blob/c440070dcb3c5bcc6ddd2227cf7adbef18573b6a/docs/phase5-readiness.md).

Until that is settled, the honest statement is: **orientation power is 0.65 and the cause is not yet
identified.** It is not cross-talk (the off-diagonal localization below shows nothing projection-associated
for orientation's own statistic), it is not latent dimensionality, and it is not a weak construction.

## Off-diagonal responses, localized

Off-diagonal rejection is *not* interpreted as estimator cross-talk by default. Each response is localized
against the checkpoint where it first becomes material, on a scale-free quantity compared only against that
same checkpoint's own zero-effect null (`report/phase4_localization.csv`).

| Mode → statistic | First material checkpoint | Reading |
|---|---|---|
| magnitude → angle | population (standardized) | construction-present |
| shape → angle | population (standardized) | construction-present |
| shape → delta | population (standardized) | construction-present |
| orientation → shape | **PLS latent** | projection-associated |
| orientation → delta | never material | not material |
| magnitude → shape | never material | not material |
| translation → anything | never material | not material |

Most off-diagonal responses are already present in the analytic population geometry, before any sampling,
preprocessing, or projection: they are properties of the constructions Phase 2 documented as mixed, and
demanding purity from them would fail a correct estimator. Magnitude's population `angle` is 39.9° yet its
`angle` rejection rate is 0.01 — the estimator is *more* specific than the construction, which is the
opposite of cross-talk.

One response is genuinely projection-associated: **orientation → shape** first becomes material at the PLS
latent checkpoint (at effects 0.75 and 1.00), where the population and observed shape values are immaterial.
Correspondingly, orientation's `shape` rejection rate is 0.97–1.00. Phase 5 should treat this as an open
question about the latent representation, not as evidence about the shape estimator. These labels describe
*where* a response appears, not what caused it.

## PLS representation

Component selection was extremely stable: the modal selected latent dimensionality was **3** in every one of
the 19 cells, with an observed range of 2–3 and mean CV AUROC 1.00 (`report/phase4_pls_selection.csv`).
Nothing in the pilot depends on a fragile component choice.

That stability has a cause worth stating: with four stages, a stage-supervised PLS-DA has at most three
informative components, and the double CV saturates there. The latent space is sized to separate stage
centroids, which is exactly what the architecture intends — but it is not sized to preserve the *group*
orientation contrast, and the attribution diagnostics quantify the cost.

## Attribution diagnostics

All 400 eligible orientation replicates produced valid diagnostics (availability 1.00), with 423 truth
drivers per replicate among 658 aligned features. Attribution cost 0.45% of total compute.

| Metric (observed component, effect 1.00) | Value |
|---|---:|
| Top-20 precision against generator truth | **1.00** |
| Top-20 recall | 0.066 |
| Bootstrap sign stability | 0.87 |
| Observed-vs-PLS-captured cosine | 0.07 |
| PLS-captured norm ratio | 0.06 |
| Top-20 precision, **PLS-captured** component | 0.15 |

Two findings stand out.

**The observed-space attribution is exact.** Every one of the top 20 observed features is a true generator
driver, at every effect size and every transition. Recall is low only because 20 identifiers are being drawn
from 423 drivers. Bootstrap sign stability of 0.80–0.87 means the reported directions are stable within a
replicate.

**The PLS reconstruction retains almost none of the orientation contrast.** The cosine between the observed
directional contrast and its PLS-captured reconstruction is 0.05–0.08, and the captured norm is ~6% of the
observed. Its top-20 precision against truth is 0.01–0.16 — near chance at small effects. Reconstructing a
~660-dimensional directional contrast through a rank-3 model cannot do better, so this is arithmetic, not a
bug; but it means the `pls_captured` view must **not** be presented as a driver list. For interpretation,
the observed component is the trustworthy one.

Cross-replicate top-20 overlap is low (Jaccard 0.02–0.05) and cross-replicate sign agreement is ~0.68. This
is expected by construction and is **not** an instability finding: each replicate index draws a fresh set of
differential indicators, so the true driver set genuinely differs between replicates. Matched seeds pair
cells *within* a replicate index, not across them. A cross-replicate stability claim would require holding
the driver set fixed, which this design deliberately does not do.

## Limitations

- 100 replicates give a standard error of about 0.05 at a rate of 0.5; small differences between adjacent
  power points are not resolvable, which is why monotonicity is judged within two combined standard errors.
- The four modes' `0.00` points come from one shared anchor cell. This is deliberate — at zero requested
  effect the generator ignores the mode, so separate cells would be byte-identical — but it means the four
  null rows are one measurement and must not be counted as independent evidence.
- Matched seeds make requested-effect comparisons paired within a replicate index. Comparisons *across*
  replicate indices remain unpaired.
- The generator-truth definition covers features whose group-stage differential mean change differs between
  groups, including propagated CpG→gene→protein effects. It does not model correlation-induced drivers, so
  precision against it is a lower bound on biological plausibility, not a causal claim.
- Attribution is descriptive. It carries no p-values and no causal interpretation, per the Phase 3 contract.
- `n_jobs` is part of each cell's parameter signature because RRPP seeds one RNG stream per worker. This run
  used the config's `n_jobs=1` throughout; parallelism came from concurrent shards. Overriding `--n-jobs`
  would change the realized permutation draws and break resume.
- The pilot executed from the working tree that became commit `073bb2f`. Record generation used exactly that
  code; the localization scaling fix and the `n_jobs` guard landed afterwards and affect only report
  generation, and the reports here were regenerated after them and verified byte-identical.

## Decision and recommended next work

**Do not launch the Phase 5 paper-grade study yet.** More Monte Carlo precision would not resolve the
orientation shortfall — it would only measure a flat 0.65 curve more accurately.

Before Phase 5:

1. **Characterize the orientation noise floor.** Translation's 16.8° observed angle at a true 0° is the
   binding constraint. Establish how the floor scales with samples per group-stage cell and with feature
   count, and pick a Phase 5 design point where orientation is adequately powered — or revise the power
   target to what the estimator can achieve at a defensible sample size.
2. **Investigate orientation → shape at the PLS checkpoint.** It is the one response that first appears at
   projection. Determine whether it is a property of a rank-3 stage-supervised space or of the shape
   statistic measured within it.
3. **Reconsider latent dimensionality for group contrasts.** Selection saturates at `n_stages - 1` because
   the supervision target is the stage label. A space sized for stage separation retains ~6% of the group
   orientation contrast. If orientation is a primary estimand, the latent space that measures it may need a
   different sizing criterion — this is a design question, not a bug.
4. **Present only the observed component in driver reports.** Its top-20 precision is 1.00; the
   PLS-captured component's is 0.15.

Magnitude, shape, and Type I calibration need no further pilot work. They are ready for Phase 5 as soon as
the orientation question is settled.

## Reproduction

From the repository root, with the committed configuration:

```bash
# 16 resumable shards (do NOT pass --n-jobs: it changes the permutation draws
# and the cell parameter signature).
for i in $(seq 0 15); do
  uv run python scripts/run_study_shard.py \
    --config examples/trajectory_power_study/phase4_pilot_100x199.json \
    --out-dir results/phase4-2026-08-27 \
    --shard-index "$i" --n-shards 16 \
    --error-policy record &
done
wait

uv run python scripts/motco_study.py merge \
  --out-dir results/phase4-2026-08-27

uv run python scripts/motco_study.py report \
  --config examples/trajectory_power_study/phase4_pilot_100x199.json \
  --out-dir results/phase4-2026-08-27
```

The per-shard and merged JSONL (61 MB) are gitignored as regenerable; the `report/` outputs and
`PROVENANCE.txt` are committed, and every claim above traces to one of them under
`results/phase4-2026-08-27/report/`:
`phase4_operating.csv`, `phase4_geometry.csv`, `phase4_pls_selection.csv`, `phase4_attribution.csv`,
`phase4_localization.csv`, `phase4_gate.csv`, and `phase4_gate_decision.json`, with figures
`phase4_geometry.png`, `phase4_selected_components.png`, and `phase4_attribution_stability.png` alongside
the existing specificity, Type I, and power outputs.
