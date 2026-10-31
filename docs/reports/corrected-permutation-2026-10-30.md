# Corrected-permutation check — 2026-10-30

> This report validates the residual-corrected RRPP null that MOTCO 0.13.0 offers as the opt-in
> `residual_correction` (openspec change `residual-df-corrected-permutation`, spec `corrected-permutation-validation`).
> It is the follow-up that the [case-study configuration diagnosis](case-study-diagnosis-2026-10-30.md) named with
> its `residual_df_shrinkage` verdict, and roadmap item 8.
>
> Every number below is read from `results/corrected-permutation-2026-10-30/`:
> - `corrected_permutation_decision.json`, `size_readings.csv`, `decide.json`;
> - per run in `type-i-astrocyte-delta0.1/`, `type-i-astrocyte-delta0.05/`, `type-i-microglia-delta0.1/`,
>   `type-i-microglia-delta0.05/` and `power/`: `corrected_permutation_rates.csv`, `merge.json`, `settings.json` and
>   `PROVENANCE.txt`, plus `corrected_permutation_pvalues.csv` for the Type I runs and
>   `corrected_permutation_contrasts.csv` for the power run;
> - the raw accounting `sacct.txt`, `finish_git.txt` and `finish_sha256.txt`;
> - the smoke records `smoke-*/` and `smoke-decision/`.
>
> `records/` and the power run's 2.6 MB `corrected_permutation_pvalues.csv` stay on `ing`; each `PROVENANCE.txt` names
> them with their sha256. Rates are rejection rates at p ≤ 0.05 over R = 500 replicates at 999 permutations.

**Run.**
- Five SLURM arrays on `ing` (partition `512x1024`, n11 excluded), from the clean worktree of
  `feat/residual-df-corrected-permutation` at `a2f1a1f`. Its git tree `4d8a22fb…` equals the workstation's, and the
  change's tests passed there.
- Jobs 930470, 930520, 930570 and 930620 ran the four Type I runs (50 shards each) and job 930670 the power run
  (140 shards). All 340 tasks ended `COMPLETED 0:0`, between 14:00Z and 14:59Z on 2026-10-30.
- 1,000 units per Type I run and 7,000 for power, each analysed three times (no correction, `df`, `leverage`). There
  were 0 arm failures, and no unit's arms disagreed on latent scores, rank or statistics.
- 164.1 elapsed core-hours (planned 160); the units' own timers sum to 123.1 h.
- The afterany job 930817 fast-forwarded the worktree to `ea9bacc`, which adds only the decision code's clustered
  sensitivity (`finish_git.txt`: nothing changed under `src/` or in the construction scripts). It then ran `merge` per
  run and `decide`. The five settings hashes in the decision file equal those of the five `settings.json` files.

**Verdict: `adopt_df`**, read verbatim from `corrected_permutation_decision.json` (`"indicative": false`).
- `"holds": {"df": true, "leverage": true}`, and `"failing"` is empty for both corrections.
- The tie-break found no significant paired difference between them on any target pair, so the rule chose `df`.
- `"uncorrected_excess_reproduced": false`. The fresh seed family did not reproduce the uncorrected `delta` excess;
  the uncorrected arm failed on pooled `shape` and on one `angle` reading instead.

**In one paragraph.**
- **The design.** The case study's configuration was rebuilt on a seed family no committed check uses:
  - six covariate columns, the rank floor of 2, net displacement and the weak stage signal;
  - both case-study cohort tables at δ = 0.1 and 0.05, with the anchor and translation controls;
  - one power run at the Phase 6 stage effect.
  
  Every unit ran `analyze_trajectories` three times on the same data, latent basis and RRPP seed, under no correction,
  `df` and `leverage`.
- **Size.** Both corrections hold every reading and every pooled rate. `df` reads 0.029–0.058 pooled per table and
  statistic against 0.0597, and at most 32 / 500 per reading against 39 / 500.
- **The uncorrected null.** It is liberal on pooled `shape` in both tables (0.0685 and 0.0675) and on the power
  design's anchor `angle` (43 / 500), but not on `delta` (0.046 and 0.045).
- **Cost.** The correction costs 0.014–0.054 of power on the target statistics. On `angle` and `shape` (0.024–0.034)
  that is of the order of the uncorrected arm's own excess at the controls.
- **`leverage`.** It is indistinguishable from `df`: at most 3 discordant units in either direction on any target.
- **Consequence.** `sead_sexdiff` re-issues both cell types' primaries with `--residual-correction df`, pinned to
  0.13.0, in its own change. MOTCO documents `df` as validated at this configuration and keeps it opt-in.

## 1. Construction

`scripts/corrected_permutation_check.py` (design D5–D7) builds every unit with the parent check's code
(`case_study_type_i_check.build_unit`) and settings, with four overrides: base seed 2031, the named cohort table, the
profile's δ, and the profile's cells.

**Fixed across runs** (each run's `settings.json`):
- the Phase 6 replicate: three stages, methylation + expression, joint magnitude and `p_dmp` 0.1;
- the confounded covariate `z = 0.5·g + 1·s + 2·g·s + ε`, injected at `s_z = 10` into a fifth of the features along a
  fixed direction;
- five i.i.d. N(0, 1) noise covariates `x1`–`x5`, so the trajectory model declares six covariate columns and p_red = 10;
- the rank floor of 2, net displacement, the harness CV knobs, 500 replicates and 999 permutations.

**Fresh family** (design D7). Base seed 2031 is used by no committed check. Every seed derives from it and the replicate
index, including the injection direction, the covariates, the generator and RRPP. The diagnosis chose the remedy on
the parent's base seed 2026, so a test on those replicates would have been circular.

**Runs.**

| run | cohort table | n | κ | δ | cells |
|---|---|---|---|---|---|
| `type-i-astrocyte-delta0.1` | `[[11, 10, 28], [9, 10, 12]]` | 80 | √(80/70) = 1.069 | 0.1 | anchor, translation 1.0 |
| `type-i-astrocyte-delta0.05` | as above | 80 | 1.069 | 0.05 | anchor, translation 1.0 |
| `type-i-microglia-delta0.1` | `[[11, 10, 24], [9, 11, 11]]` | 76 | √(76/66) = 1.073 | 0.1 | anchor, translation 1.0 |
| `type-i-microglia-delta0.05` | as above | 76 | 1.073 | 0.05 | anchor, translation 1.0 |
| `power` | astrocyte | 80 | 1.069 | 2.0 | anchor, translation 1.0, magnitude 0.02–0.25, orientation and shape 0.25–1.0 (14 cells) |

The two tables are the `cohort_30` tables of the case study's `result.json` files (sha256 recorded in each
`settings.json`). The two controls of a replicate, and the power run's 14 cells, share that replicate's seeds.

**Arms.** Each unit calls `analyze_trajectories` with identical parameters except `residual_correction` ∈
{`None`, `"df"`, `"leverage"`}. The unit then asserts that the three latent-score arrays, the selected ranks and the
statistics are exactly equal. A difference would raise and fail the shard; none occurred.

**Calibration** (descriptive, from the uncorrected arm at the anchor; `calibration` in the decision file):

| run | mean held-out CV AUROC (band 0.60–0.75) | median \|ρ\| with `z` (≥ 0.8) | selected rank 2 / 3 / ≥ 4 |
|---|---|---|---|
| astrocyte, δ = 0.1 | 0.708 | 0.929 | 341 / 150 / 9 |
| astrocyte, δ = 0.05 | 0.693 | 0.929 | 386 / 111 / 3 |
| microglia, δ = 0.1 | 0.705 | 0.937 | 322 / 163 / 15 |
| microglia, δ = 0.05 | 0.689 | 0.936 | 363 / 132 / 5 |
| power | 1.000 (outside the band, as expected at δ = 2.0) | 0.802 | 290 / 210 / 0 |

The construction was calibrated on the astrocyte table and not re-tuned. The microglia table also reads inside the
band.

**Rule** (design D8; `rule` in the decision file). A correction holds size when both tiers hold:
- **pooled**, per cohort table × statistic over the table's four Type I readings: rate ≤ α + 2·√(α(1 − α)/2000) =
  0.0597, i.e. at most 119 / 2,000;
- **per reading**, over the four Type I runs' two controls and the power run's two controls, each for three statistics
  (30 readings): rate ≤ α + 3·√(α(1 − α)/500) = 0.0792, i.e. at most 39 / 500.

The uncorrected arm goes through the same rule descriptively. A replicate-clustered pooled sensitivity was added on
2026-10-30, before any paper-grade result. It is descriptive and never changes the verdict (design, "Smoke and shard
plan").

**Smoke** (`smoke-*`, R = 20 × 99 permutations). Every unit completed with 0 failures, and the indicative verdict was
`no_correction_holds`. The only failing entries were the pooled `shape` readings, 9–10 / 80, in every arm including
the uncorrected one; at R = 20 that was read as noise. At R = 500 the corrected arms hold `shape` and the uncorrected
arm does not (§2).

## 2. Size

### Per reading

Rejections out of 500 as none / `df` / `leverage`, from `size_readings.csv` (tier `reading`). **Bold** marks a reading
above 39 / 500.

| run | control | `delta` | `angle` | `shape` |
|---|---|---|---|---|
| astrocyte, δ = 0.1 | anchor | 19 / 13 / 13 | 33 / 24 / 23 | 33 / 30 / 30 |
| astrocyte, δ = 0.1 | translation | 26 / 18 / 18 | 24 / 20 / 20 | 35 / 29 / 29 |
| astrocyte, δ = 0.05 | anchor | 21 / 12 / 11 | 30 / 23 / 24 | 36 / 31 / 32 |
| astrocyte, δ = 0.05 | translation | 26 / 19 / 20 | 25 / 18 / 21 | 33 / 26 / 26 |
| microglia, δ = 0.1 | anchor | 20 / 14 / 14 | 27 / 23 / 24 | 30 / 22 / 22 |
| microglia, δ = 0.1 | translation | 25 / 16 / 17 | 30 / 27 / 29 | 37 / 28 / 27 |
| microglia, δ = 0.05 | anchor | 20 / 12 / 13 | 27 / 22 / 22 | 33 / 27 / 27 |
| microglia, δ = 0.05 | translation | 25 / 16 / 17 | 30 / 28 / 27 | 35 / 27 / 26 |
| power (astrocyte, δ = 2.0) | anchor | 23 / 17 / 17 | **43** / 32 / 34 | 32 / 23 / 22 |
| power (astrocyte, δ = 2.0) | translation | 28 / 20 / 21 | 35 / 29 / 30 | 35 / 28 / 28 |

- `df`'s largest reading is 32 / 500 (power anchor `angle`), and `leverage`'s is 34 / 500 at the same reading. Both
  are within the parent check's tighter k = 2 bound of 34 / 500 as well.
- The uncorrected arm exceeds once: power anchor `angle`, 43 / 500 = 0.086.

### Pooled per cohort table

Rejections out of 2,000 (rate), from `size_readings.csv` (tier `pooled`). **Bold** marks a rate above 0.0597.

| table | statistic | none | `df` | `leverage` |
|---|---|---|---|---|
| astrocyte | `delta` | 92 (0.046) | 62 (0.031) | 62 (0.031) |
| astrocyte | `angle` | 112 (0.056) | 85 (0.0425) | 88 (0.044) |
| astrocyte | `shape` | **137 (0.0685)** | 116 (0.058) | 117 (0.0585) |
| microglia | `delta` | 90 (0.045) | 58 (0.029) | 61 (0.0305) |
| microglia | `angle` | 114 (0.057) | 100 (0.05) | 102 (0.051) |
| microglia | `shape` | **135 (0.0675)** | 104 (0.052) | 102 (0.051) |

### Replicate-clustered sensitivity

Rate / bound (design effect), from `size_readings.csv` (tier `pooled_clustered`). The bound is
α + 2·√(α(1 − α)/R_eff), where R_eff is 2,000 divided by the design effect of the rate over 500 replicate clusters.

| table | statistic | none | `df` | `leverage` |
|---|---|---|---|---|
| astrocyte | `delta` | 0.046 / 0.0632 (1.84) | 0.031 / 0.0630 (1.77) | 0.031 / 0.0630 (1.77) |
| astrocyte | `angle` | 0.056 / 0.0632 (1.83) | 0.0425 / 0.0631 (1.80) | 0.044 / 0.0632 (1.84) |
| astrocyte | `shape` | **0.0685 / 0.0629 (1.75)** | 0.058 / 0.0630 (1.79) | 0.0585 / 0.0631 (1.80) |
| microglia | `delta` | 0.045 / 0.0633 (1.86) | 0.029 / 0.0632 (1.84) | 0.0305 / 0.0633 (1.86) |
| microglia | `angle` | 0.057 / 0.0638 (1.99) | 0.05 / 0.0638 (2.00) | 0.051 / 0.0639 (2.04) |
| microglia | `shape` | **0.0675 / 0.0638 (1.99)** | 0.052 / 0.0637 (1.98) | 0.051 / 0.0636 (1.94) |

- The design effects are 1.75–2.04, so the four readings of a table carry about 980–1,140 independent replicates'
  worth of information, not 2,000. The smoke had predicted that (design, "Smoke and shard plan").
- The verdict does not depend on it. `"clustered_sensitivity_holds"` is `{"none": false, "df": true,
  "leverage": true}`, the same split as the predeclared rule.

## 3. The uncorrected null

**`uncorrected_excess_reproduced: false`.** The flag asks whether the uncorrected pooled `delta` exceeds 0.0597 in
either table. It reads 0.046 (astrocyte) and 0.045 (microglia). On the parent's family the same model read 0.066
pooled (`z_noise5` in the [diagnosis](case-study-diagnosis-2026-10-30.md), §2).

As predeclared (design, Risks), the correction's size verdict stands either way, and no rule was re-tuned. The
evidence that column count drives the `delta` excess remains the diagnosis's within-family contrast C2 (p = 0.015).
This check runs every arm at the same six covariates, so it has no column contrast and cannot re-test that mechanism.

**The uncorrected null is still liberal at this configuration, on other statistics.** The descriptive `uncorrected`
block lists three failing entries:
- pooled `shape` 0.0685 (astrocyte) and 0.0675 (microglia) against 0.0597. They also fail the clustered sensitivity.
  The parent family read 0.0615 for the same model (diagnosis, §2).
- power anchor `angle`, 43 / 500 = 0.086 against 39 / 500, at the Phase 6 stage effect.

Every corrected arm holds all three: `df` reads `shape` 0.058 and 0.052 pooled and power anchor `angle` 32 / 500. So the
narrowed Freedman–Lane null shows up on different statistics in different seed families and designs:
- `delta` on the parent family;
- `shape` on this family's Type I runs;
- `angle` at this family's strong-signal anchor.

The correction removes the excess in each case, because it widens the null for every statistic at once.

## 4. Power cost

From `power/corrected_permutation_rates.csv`: rejection rates on each mode's target statistic, as none / `df` /
`leverage`, by effect.

| mode / statistic | effect 1 | effect 2 | effect 3 | effect 4 |
|---|---|---|---|---|
| magnitude / `delta` | 0.02: 0.082 / 0.056 / 0.056 | 0.05: 0.238 / 0.19 / 0.192 | 0.1: 0.536 / 0.482 / 0.478 | 0.25: 0.94 / 0.926 / 0.928 |
| orientation / `angle` | 0.25: 0.182 / 0.154 / 0.154 | 0.5: 0.29 / 0.256 / 0.26 | 0.75: 0.436 / 0.412 / 0.412 | 1.0: 0.732 / 0.706 / 0.706 |
| shape / `shape` | 0.25: 0.1 / 0.074 / 0.074 | 0.5: 0.148 / 0.118 / 0.116 | 0.75: 0.208 / 0.182 / 0.182 | 1.0: 0.296 / 0.264 / 0.262 |

Pooled over each mode's four effects (2,000 paired readings in 500 replicate clusters;
`power/corrected_permutation_contrasts.csv`, scope `mode`, role `descriptive`). "None only" counts the readings the
uncorrected arm rejects and the correction does not, and p is the clustered one-sided sign-flip p-value for that
direction. Targets are in bold.

| mode | statistic | none | `df` | none only | p | `leverage` | none only | p |
|---|---|---|---|---|---|---|---|---|
| magnitude | **`delta`** | 0.449 | 0.4135 | 71 | 8.5e-22 | 0.4135 | 71 | 8.5e-22 |
| magnitude | `angle` | 0.0625 | 0.046 | 33 | 4.8e-07 | 0.0465 | 32 | 9.5e-07 |
| magnitude | `shape` | 0.044 | 0.0285 | 31 | 1.9e-06 | 0.0285 | 31 | 3.8e-06 |
| orientation | `delta` | 0.2965 | 0.256 | 81 | 8.3e-25 | 0.2575 | 78 | 6.6e-24 |
| orientation | **`angle`** | 0.41 | 0.382 | 56 | 4.4e-16 | 0.383 | 54 | 8.9e-16 |
| orientation | `shape` | 0.4725 | 0.441 | 63 | 2.2e-19 | 0.4425 | 60 | 1.7e-18 |
| shape | `delta` | 0.273 | 0.232 | 82 | 3.3e-24 | 0.231 | 84 | 8.3e-25 |
| shape | `angle` | 0.056 | 0.0435 | 25 | 3.0e-08 | 0.0435 | 25 | 3.0e-08 |
| shape | **`shape`** | 0.188 | 0.1595 | 57 | 5.6e-17 | 0.1585 | 59 | 1.4e-17 |

- **The cost is one-directional.** No reading in any cell is rejected by a correction and not by the uncorrected arm
  (every `df`-only and `leverage`-only count in the contrast file is 0). That is expected: the observed statistics
  are identical, and only the null widens.
- **`angle` and `shape`.**
  - The orientation `angle` power falls by 0.024–0.034 per effect, from 0.732 to 0.706 at e = 1.0.
  - The shape `shape` power falls by 0.026–0.032, from 0.296 to 0.264 at e = 1.0.
  - Off-target `angle` and `shape` fall too, towards or below α: magnitude `angle` 0.0625 → 0.046 and magnitude
    `shape` 0.044 → 0.0285.
- **`delta`.** The magnitude `delta` power falls by 0.026, 0.048, 0.054 and 0.014 at e = 0.02, 0.05, 0.1 and 0.25. That
  is the largest absolute cost, in the rising part of the power curve.
- **The comparison arm is itself liberal here.** At the power run's controls the uncorrected arm rejects more than
  `df` by 0.012 / 0.022 / 0.018 (anchor `delta` / `angle` / `shape`) and 0.016 / 0.012 / 0.014 (translation). The
  `angle` and `shape` costs above are of that order. A rate difference against an arm whose null is too narrow
  therefore overstates the power a size-valid test gives up. It is reported as measured, against no correction, as
  the spec requires.

The diagnosis saw the same direction at the null on the parent family: `z_noise5_corrected` read `angle` 0.042 and
`shape` 0.0465, against 0.0455 and 0.0615 uncorrected.

## 5. Tie-break

`leverage` against `df` on the power run's target pairs, pooled over each mode's four effects (decision file
`tie_break`; contrast file role `tie_break`):

| target | `leverage` rejections | `df` rejections | `leverage` only | `df` only | p (gain) | p (loss) |
|---|---|---|---|---|---|---|
| magnitude / `delta` | 827 | 827 | 3 | 3 | 0.656 | 0.656 |
| orientation / `angle` | 766 | 764 | 2 | 0 | 0.25 | 1.0 |
| shape / `shape` | 317 | 319 | 1 | 3 | 0.938 | 0.313 |

No gain or loss reaches p ≤ 0.05, so `"choice": "df"` (design D9). The two corrections disagree on at most 3 of 2,000
readings in either direction on any target, and their size readings differ by at most 3 / 500 per reading and 3 / 2,000
pooled.

The per-sample multipliers do differ from κ; each analysis's range is in its record's
`measurement.residual_correction`, kept in `records/`. Even so, at these tables the two corrections almost never reach
different rejection decisions.

## 6. Verdict and consequences

Read verbatim from `corrected_permutation_decision.json`:

```json
{"verdict": "adopt_df", "indicative": false, "holds": {"df": true, "leverage": true},
 "failing": {"df": [], "leverage": []}, "uncorrected_excess_reproduced": false}
```

**Consequence** (the decision file's `consequence`): "sead_sexdiff re-issues both cell types' primaries with
--residual-correction df, pinned to the release that contains it (0.13.0), in its own change; MOTCO documents df as
validated at this configuration and keeps it opt-in".

**For `sead_sexdiff`.**
- The covariate-adjusted primary, withheld since the [case-study Type I check](case-study-type-i-2026-10-29.md), is
  re-issued for both cell types with `--residual-correction df` under MOTCO 0.13.0, in its own change.
- Under `df` its `delta`, `angle` and `shape` p-values are size-controlled at this configuration: pooled 0.029–0.058
  per table and statistic, at most 32 / 500 per reading.
- `delta` is conservative under `df` here (0.029 and 0.031 pooled). A `delta` non-rejection therefore stays weak
  evidence, as before.
- The earlier caveats stand: the `angle` null is too wide to inform a non-rejection, and the CV-selected rank is
  reported beside any `angle` result.
- The `stage` latent-adjusted arm is not touched by this check and stays inadmissible at this configuration
  (diagnosis, `residualized_basis`).

**For MOTCO.**
- `residual_correction="df"` (`--residual-correction df`) is documented as validated at the case study's
  configuration. It stays opt-in, and without it every output equals 0.12.0's.
- `leverage` holds size too, and it is kept. It is not the validated choice, because it earned no paired advantage.
- Whether `df` should become the analysis default is a separate, later change. A default flip changes every
  covariate-adjusted output and needs its own fixture re-issue (tasks, "Workflow follow-up").

## 7. Limitations and scope

- **`cohort_50` tables not run.** The case study's sensitivity tables, for example microglia `cohort_50` at n = 63
  with κ = 1.090, were not run. The shrinkage grows with p_red / n, so the uncorrected null is expected to be more
  liberal there and the correction to matter more. The correction's size there is inferred, not measured.
- **Additive-linear covariates.** `z` enters the features additively along one fixed direction, and the trajectory
  model linearly. The five noise covariates are i.i.d. normals, unlike the case study's age, education, PMI, RIN and
  library sizes. The correction depends on n and the reduced rank only, but its size at the real covariates is
  inferred.
- **Structural-analogue blocks.** The blocks are InterSIM-calibrated methylation and expression analogues
  (498 measured features), not the case study's TMM log2-CPM blocks of 5,000 features each.
- **Calibration on the astrocyte table only.** δ and `s_z` were calibrated on the astrocyte table by the parent check
  and not re-tuned. The microglia runs read inside the band (§1), so this did not bite, but it was not designed for.
- **Power on one table and one strength.** The power run uses the astrocyte table at δ = 2.0, where the statistics
  have power. The weak-signal runs hold only null controls, so the cost at the case study's stage signal is not
  measured.
- **The `delta` excess was not replicated.** The mechanism rests on the diagnosis's C2 (§3). The uncorrected arm's
  pooled `delta` on this family is at α, which means that `df` is conservative for `delta` here.
- **What stays on `ing`.** `records/` (4.0 MB per Type I run, 26 MB for power) and the power run's 63,000-row
  `corrected_permutation_pvalues.csv` (2.6 MB) stay there, under
  `/home1/tgonza/MOTCO/results/corrected-permutation-2026-10-30/` with the SLURM logs. They were moved there from the
  run's worktree on 2026-10-31, and their digests were rechecked. Each `PROVENANCE.txt` names them with their sha256. `corrected_permutation_check.py merge` regenerates every committed table from `records/`.
