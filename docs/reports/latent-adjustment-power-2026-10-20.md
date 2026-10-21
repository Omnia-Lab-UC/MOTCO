# Latent-adjustment power comparison — 2026-10-20

> Paper-grade power comparison of covariate-aware latent construction against the 0.8.0 covariate-blind basis at
> the Phase 6 SEA-AD design (`openspec/changes/covariate-aware-latent-construction/`, design D9, spec
> `latent-adjustment-validation`). It is the companion of the [Type I check](latent-adjustment-type-i-2026-10-20.md),
> which admitted `stage` conditioning. Every number below is read from `results/latent-adjustment-power-2026-10-20/`
> (`latent_adjustment_power_decision.json`, `latent_adjustment_power_rates.csv`,
> `latent_adjustment_power_pvalues.csv`, `merged.jsonl`). Rates are rejection rates at p ≤ 0.05 over the available
> replicates of 500, at 999 permutations. The Monte Carlo SE is at most 0.024.

**Run:** 14,000 units of three analyses on `ing` (SLURM array 923444, 200 shards, 158 recorded core-hours, no
failed task), from `979657c` plus this change's working tree; file hashes are in `PROVENANCE.txt`. **Advisory
verdicts:** `stage` `candidate_fails`, `group_stage` `candidate_fails`.

**In one paragraph.** `stage` conditioning gains a lot of orientation power. At the top effect, `angle` rejects
0.922 against the covariate-blind basis's 0.784 with a covariate-aligned nuisance axis, and 0.931 against 0.642
without one. It fails the predeclared rule for three reasons:

- it loses `shape` power at the top shape effect when the nuisance axis is present (0.243 against 0.428);
- it loses power for small magnitude effects and on `delta` under shape effects, without the nuisance axis
  (for example 0.400 against 0.520 at magnitude e = 0.1);
- its translation control rejects `angle` at 0.085–0.092, above the bound. The covariate-blind basis is above its
  own bound there too (0.076–0.086).

`group_stage` conditioning loses most of the orientation and shape power in both columns: `angle` 0.344 and `shape`
0.054 at the top effects with the nuisance axis. With the Type I check's `angle` excess, it is ruled out. No arm in
either column ever selected rank 1. Latent adjustment is therefore not a uniform power improvement. It trades shape
and small-magnitude sensitivity for orientation sensitivity in this construction, which §6 turns into conditions for
the case study.

## 1. Construction

- **Cells.** The zero-effect anchor, the translation control at e = 1.0, and three modes at the Phase 6 effect
  levels:
  - magnitude (joint) at e = 0.02, 0.05, 0.10 and 0.25, below its saturation;
  - orientation at e = 0.25, 0.5, 0.75 and 1.0;
  - shape at e = 0.25, 0.5, 0.75 and 1.0.
  - The design, blocks and estimator are the [Type I check](latent-adjustment-type-i-2026-10-20.md)'s (§1).
- **Columns.**
  - `nuisance`: the Type I check's calibrated injection (`s_z = 40`, `s_c = 2`).
  - `clean`: `s_z = s_c = 0`. The same eight covariates are declared, but none affects the features.
- **Matched seeds.** The generator, covariate and RRPP seeds depend on the replicate index only. Within a
  replicate, the two columns share the generated blocks and covariates and differ only in the injection. Every cell
  starts from the same baseline trajectory, and the arms share the dataset and the RRPP seed. The nuisance column's
  anchor is therefore the Type I check itself, and its rates reproduce that report's to the third decimal.
- **Arms.** `baseline` (covariates in the trajectory model only, covariate-blind basis), `stage` and `group_stage`.
- **Rule** (fixed in the hashed settings before the run; the study engine's `column_decision` convention). Each
  residualized arm is compared against `baseline`. The pooled SE is `√(SE_arm² + SE_baseline²)` and the bound is
  0.05 + 2·√(0.05·0.95/R).
  - **gain:** `nuisance` column, orientation/`angle` at e = 1.0: `rate_arm − rate_baseline > 2 · pooled SE`;
  - **protected:** `nuisance` column, magnitude/`delta` at e = 0.25 and shape/`shape` at e = 1.0:
    `rate_baseline − rate_arm ≤ 2 · pooled SE`;
  - **no loss:** `clean` column, every statistic at every effect of every mode: the same loss tolerance;
  - **controls:** both columns, anchor and translation, every statistic at or below the bound.

  D9 also listed magnitude/`delta` and shape/`shape` as gain targets. The spec's single gain target is the one
  applied, because magnitude/`delta` is at 0.92 or above in every arm and column at its top effect and has little room to gain.

The same 52 (`stage`) and 38 (`group_stage`) replicates are unavailable in every cell. Each loses a factor level
from a CV training fold, which depends only on the covariates and the folds (Type I report §2). So R = 448 for
`stage`, 462 for `group_stage` and 500 for `baseline` throughout.

## 2. Power on the target statistic

From `latent_adjustment_power_rates.csv`, rate for the mode's own statistic (SE ≤ 0.024).

| mode | e | nuisance: baseline | nuisance: `stage` | nuisance: `group_stage` | clean: baseline | clean: `stage` | clean: `group_stage` |
|---|---|---|---|---|---|---|---|
| magnitude/`delta` | 0.02 | 0.072 | 0.042 | 0.093 | 0.076 | 0.042 | 0.084 |
| | 0.05 | 0.226 | 0.121 | 0.229 | 0.180 | 0.121 | 0.240 |
| | 0.10 | 0.560 | 0.417 | 0.591 | 0.520 | 0.400 | 0.597 |
| | 0.25 | 0.940 | 0.969 | 0.976 | 0.924 | 0.960 | 0.972 |
| orientation/`angle` | 0.25 | 0.208 | 0.138 | 0.167 | 0.128 | 0.156 | 0.165 |
| | 0.50 | 0.418 | 0.263 | 0.229 | 0.170 | 0.254 | 0.221 |
| | 0.75 | 0.584 | 0.549 | 0.255 | 0.200 | 0.500 | 0.251 |
| | 1.00 | 0.784 | **0.922** | 0.344 | 0.642 | **0.931** | 0.325 |
| shape/`shape` | 0.25 | 0.104 | 0.078 | 0.056 | 0.066 | 0.071 | 0.054 |
| | 0.50 | 0.200 | 0.080 | 0.067 | 0.090 | 0.074 | 0.063 |
| | 0.75 | 0.332 | 0.167 | 0.056 | 0.140 | 0.181 | 0.058 |
| | 1.00 | 0.428 | 0.243 | 0.054 | 0.128 | 0.250 | 0.045 |

Controls, rate per statistic (`delta` / `angle` / `shape`). The bound is 0.0695 at R = 500, 0.0706 at 448 and
0.0703 at 462.

| control | nuisance: baseline | nuisance: `stage` | nuisance: `group_stage` | clean: baseline | clean: `stage` | clean: `group_stage` |
|---|---|---|---|---|---|---|
| anchor | 0.048 / 0.062 / 0.078 | 0.031 / 0.067 / 0.058 | 0.048 / 0.082 / 0.069 | 0.044 / 0.076 / 0.062 | 0.038 / 0.067 / 0.051 | 0.052 / 0.082 / 0.074 |
| translation | 0.052 / 0.076 / 0.070 | 0.025 / 0.092 / 0.054 | 0.065 / 0.110 / 0.106 | 0.050 / 0.086 / 0.050 | 0.025 / 0.085 / 0.060 | 0.063 / 0.104 / 0.093 |

As a cross-check, the clean column's baseline orientation/`angle` at e = 1.0 (0.642) matches the Phase 6 small-n
study's net-displacement rate at the same effect (0.648; [power comparison](net-displacement-power-2026-10-14.md)),
even though eight covariates are now declared.

## 3. The verdicts

From `latent_adjustment_power_decision.json`.

**`stage` — `candidate_fails`.**

- **gain:** passes. 0.922 against 0.784, difference 0.138 > 0.045.
- **protected:** magnitude/`delta` passes (0.969 against 0.940). shape/`shape` fails: 0.243 against 0.428, a loss
  of 0.185 against a tolerance of 0.060.
- **no loss** (clean column) fails in eight reads:
  - magnitude/`delta` at e = 0.02, 0.05 and 0.10: 0.042, 0.121 and 0.400 against 0.076, 0.180 and 0.520;
  - `delta` under orientation at e = 0.25: 0.029 against 0.066;
  - `delta` under shape at all four effects: for example 0.214 against 0.298 at e = 1.0.
  - Every other clean-column read holds. That includes shape/`shape`, which gains: 0.250 against 0.128 at
    e = 1.0.
- **controls:** translation/`angle` exceeds in both columns (0.092 and 0.085 against 0.0706). Every other control
  holds.

**`group_stage` — `candidate_fails`.** The gain fails (0.344 against 0.784) and shape/`shape` fails as a protected
pair (0.054 against 0.428). Nine clean-column no-loss reads fail, all under orientation or shape effects. Seven
controls exceed the bound: anchor `angle` in both columns at 0.082, anchor `shape` in the clean column, and
translation `angle` and `shape` in both columns.

## 4. Selected rank

Pooled over the 14 cells of each column, as a share of the available units:

| column | arm | rank 2 | rank 3 | rank ≥ 4 | rank 1 |
|---|---|---|---|---|---|
| nuisance | baseline | 0.381 | 0.513 | 0.106 | 0 |
| nuisance | `stage` | 0.048 | 0.842 | 0.110 | 0 |
| nuisance | `group_stage` | 0.010 | 0.984 | 0.006 | 0 |
| clean | baseline | 0.341 | 0.621 | 0.038 | 0 |
| clean | `stage` | 0.045 | 0.853 | 0.102 | 0 |
| clean | `group_stage` | 0.008 | 0.987 | 0.005 | 0 |

Per-cell distributions are in the `rank_distribution` column of the rates table. Residualization moves the CV
almost entirely off rank 2. **No unit in any arm or column selected rank 1.** This construction therefore does not
reproduce the SEA-AD cohorts' rank-1 selection and is not evidence about it. The separate rank-1 `angle` change
needs a construction with a single dominant axis, or real-data evidence (roadmap, Next changes item 4).

## 5. Reading

- **Residualization removes signal even when no covariate has an effect.** On six clean-column orientation
  replicates at e = 1.0 (a diagnostic run on the workstation, not part of the records), the subtracted term `C γ`
  carries about 19 % of each pooled-standardized feature's variance under `stage` and about 12 % under
  `group_stage`.
  - Under `stage`, `z` differs by group within each stage by construction (`a_g + a_gs·s`), so `γ` takes in part of
    the group difference. About a quarter of the features' between-cell (group × stage) variance is removed from
    the features the basis is fitted on.
  - Under `group_stage`, 12 % is what nine covariates fit by chance on 74 within-cell degrees of freedom. About
    40 % of it lands between cells, because `z`'s cell means run from 0.56 to 6.94 against a within-cell SD of
    about 0.9, so the estimation error of `γ` is extrapolated into a spurious cell-mean pattern.
  - The basis is then fitted, at the pooled scale, on features that carry that spurious pattern. This is consistent
    with, though not proof of, `group_stage`'s loss of orientation and shape power and its control excesses.
- **`stage`'s profile is a trade.** It finds the orientation surgery far more often in both columns, and the shape
  surgery more often in the clean column. It finds small magnitude effects and the shape surgery less often when
  the nuisance axis is present. `delta` rejects less often under shape and orientation effects. Because the
  orientation gain also appears with no nuisance effect at all, it is not simply the nuisance axis being removed.
  The basis changes whenever the covariates are confounded with group within stage, as they are here by
  construction.
- **The nuisance axis also inflates the covariate-blind basis's power.** Its `shape` rate at e = 1.0 is 0.428 with
  the axis and 0.128 without, and its orientation `angle` rate is 0.784 against 0.642. Part of `stage`'s
  protected-pair `shape` "loss" is therefore a loss against a nuisance-inflated reference, consistent with the
  baseline's anchor `shape` of 0.078 above its bound in the Type I check. The rule does not distinguish the two, and
  this report does not overrule it.
- **The translation control.** `angle` under a pure location shift is at or above the bound for every arm here,
  including the covariate-blind basis (0.076 and 0.086). The [net-displacement power comparison](net-displacement-power-2026-10-14.md)
  already placed translation near the bound (0.058–0.068) without covariates. With eight covariates this
  construction pushes all arms above it, `stage` by a further 0.016 in the nuisance column.

## 6. Scope, limits and consequence

- **Mediators.** Residualizing a variable that mediates the stage effect removes real signal. No simulation can
  tell a confounder from a mediator; this comparison injects a pure confounder. The case study must choose and
  justify its covariate set in its own proposal.
- **Extreme confounding.** The covariates here are confounded with group × stage far beyond what technical
  variables usually are: `z`'s cell means span about 7 within-cell SDs. The power costs above, and the
  `group_stage` failure, grow with that confounding. Real covariates that are less tied to the design will remove
  less signal, but this comparison does not measure how much less.
- **One design.** n = 80, two blocks, 498 features, eight covariates, one injected direction per covariate.
- **Advisory.** The verdicts are advisory; the Type I decision (`adopt_stage`) is the gate.

**Consequence for the case study.** `--latent-adjustment stage` remains size-valid (Type I report §5). It is not a
power upgrade to adopt instead of the covariate-adjusted 0.8.0 analysis:

- report it as a sensitivity analysis beside the 0.8.0 analysis on the same covariates, not as a replacement;
- read an `angle` rejection that appears only under latent adjustment in light of the translation control above,
  which is at or above the bound for both bases;
- treat a `delta` or `shape` non-rejection under latent adjustment as weaker than under the covariate-blind basis;
- report the record's `latent_adjustment` block, in particular `feature_r2`, which shows how much variance the
  residualization removed;
- do not use `group_stage`.
