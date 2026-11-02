# Operating characteristics at the ADNI design — 2026-11-02

> This report measures size and power of the real-data analysis at the ADNI case study's design (openspec change
> `adni-design-operating-characteristics`, spec `adni-design-validation`; roadmap Next changes item 9). The case study
> (`adni_sexdiff`) reads its `angle` and `shape` non-rejections as uninformative under its reading rule 3, because no
> power statement existed at its design. This check supplies one.
>
> Every number below is read from `results/adni-design-2026-11-02/`:
> - `run/adni_design_decision.json`, `run/size_readings.csv`, `run/adni_design_rates.csv`,
>   `run/adni_design_rank_strata.csv`, `run/merge.json`, `run/decide.json`, `run/settings.json`, `run/PROVENANCE.txt`;
> - `calibration/calibration.json` (every calibration iteration);
> - the smoke `smoke/` (R = 20 × 99 permutations, indicative).
>
> `run/records/` and the 1.66 MB `run/adni_design_replicates.csv` stay on `ing`; `run/PROVENANCE.txt` names them with
> their sha256. Rates are rejection rates at p ≤ 0.05 over R = 500 replicates at 999 permutations.

**Run.**
- One SLURM array on `ing` (job 932122, partition `512x1024`, n11 excluded, 175 shards of 40 units), from the clean
  worktree of `feat/adni-design-operating-characteristics` at `db2e30a`. Its git tree `3ec28936…` equals the
  workstation's, and the change's tests passed there.
- All 175 tasks ended `COMPLETED 0:0` between 15:23Z and 16:48Z on 2026-11-02. 7,000 of 7,000 units, 0 failures.
- 174.6 elapsed core-hours (planned about 200, range 105–265); the units' own timers sum to 152.7 h.
- The afterany job 932123 ran `merge` and `decide`. The decision file's settings hash `aa7964ab…` equals
  `settings.json`'s.

**Verdict: `non_rejections_uninformative`**, read verbatim from `adni_design_decision.json` (`"indicative": false`).
- `"size"`: `"holds": true`, `"exceeding": []`. All six readings sit at or below the 0.0695 bound.
- `"informative_from"` is null for `delta`, `angle` and `shape`: no target statistic reaches the 0.80 floor at any
  effect.

**In one paragraph.**
- **The design.** ADNI's sex × stage table (n = 1,501) was rebuilt with two blocks and six covariate columns of ADNI's
  kinds: a confounded numeric, a factor in two drop-first columns, a phase-like binary nuisance injected into one block,
  and two noise columns. The stage signal was calibrated to ADNI's held-out AUROC, and the run used a seed family no
  committed check uses. Every unit ran the case study's primary as issued: `analyze_trajectories` with the six
  covariates, the rank floor of 2 and net displacement.
- **Size.** It holds at both controls. The six readings are 0.034–0.068 against 0.0695.
- **Power.** None of the three statistics has useful power. Magnitude/`delta` reaches 0.530 at a doubled stage effect
  in one group (e = 1.0). Orientation/`angle` and shape/`shape` stay at the null rate (0.034–0.076) at every effect.
- **Rank.** The CV-selected rank sits at the cap, 19, in 95–100 % of units at every cell, against ADNI's 8.
- **Consequence.** `adni_sexdiff` keeps its rule 3, now citing these measured rates, in its own change. The Phase 6
  exit review (roadmap item 10) reads ADNI as a design without power for any of the three statistics.

## 1. Construction

`scripts/adni_design_check.py` composes the case-study construction by import (`case_study_type_i_check`, built on
the D10 construction of `covariate_type_i_check`):

- **Generator.** The semi-synthetic generator at ADNI's table `[[237, 310, 125], [223, 448, 158]]`, three stages,
  methylation (367 CpGs) and expression (131 genes), joint magnitude, `p_dmp` 0.1, and one weakened per-omic stage
  effect δ.
- **Covariates, in model order** (six design columns, as `covariates.design_columns` records them):
  - `z`: D10's numeric covariate, confounded with group × stage and injected along a fixed direction over 20 % of the
    joint features at `s_z` (ADNI's `age`).
  - `category[c1]`, `category[c2]`: D10's three-level covariate on the same linear predictor, declared as a factor and
    not injected (ADNI's `apoe4[1]`, `apoe4[2]`).
  - `phase`: binary, drawn per sex × stage cell at ADNI's ADNI2GO shares (F 0.603 / 0.623 / 0.424, M 0.561 / 0.507 /
    0.481; `adni_sexdiff` `results/motco_input/`). It is injected into the **expression block only** as a dense offset
    `s_p · sd_f · d_f · phase`, with `d_f` i.i.d. N(0, 1) drawn once from the base seed. It enters the model as a
    numeric 0/1 column (ADNI's `phase[ADNI2GO]`). The expression block's 131 features stand for p180's 133.
  - `x1`, `x2`: i.i.d. N(0, 1), never injected (ADNI's `education`, `statin[1]`).
- **Seeds.** Every seed derives from base seed 2032 and the replicate index, so the 14 cells of a replicate share the
  baseline draw, the covariates, the phase draw and the permutation draws, and differ only in the surgery.
- **Analysis.** One `adjusted` analysis per unit: the six covariates in the trajectory model, `min_components = 2`,
  net displacement, the harness CV knobs (3 × 4 folds, 5 repeats, `max_components` 20), 999 permutations and the
  replicate's RRPP seed. There is no latent adjustment and no residual correction, as in ADNI's primary (κ ≈ 1.003 at
  n = 1,501 would be immaterial).
- **Cells.** The anchor; `translation` at e = 1.0; magnitude at e = 0.1, 0.25, 0.5, 1.0; orientation and shape at
  e = 0.25, 0.5, 0.75, 1.0.

### Calibration

Twelve iterations on anchor units (R = 20 × 99 permutations, two arms), all in `calibration/calibration.json`, met the
three predeclared criteria at **δ = 0.02, `s_z` = 0.1, `s_p` = 1.75**. The iterations and why this one was chosen are
tabulated in the change's `design.md` (Calibration). In short:

- At D10's `s_z` = 3 the first component followed `z` (|ρ| 0.94) and the AUROC sat at 0.77–0.88.
- At `s_z` = 0.3 the AUROC kept a floor near 0.645 that δ barely moved.
- At `s_z` = 0.1 the floor fell, and the phase offset steered the basis without adding stage signal: raising `s_p`
  from 1.5 to 2.5 moved the AUROC from 0.602 to 0.604.

The calibrated construction, as realized at the paper-grade run's 500 anchor units (`calibration_realized` in the
decision), beside ADNI's primary:

| quantity | construction (anchor, R = 500) | ADNI primary |
|---|---|---|
| held-out CV stage AUROC, mean | 0.606 (q05–q95 0.576–0.633; band 0.55–0.65) | 0.600 (0.580–0.616 across CV-selected analyses) |
| first latent component \|Spearman\| with phase, median | 0.578 (threshold 0.4) | 0.542 |
| first latent component \|Spearman\| with the confounded numeric, median | 0.078 | 0.061 (`age`) |
| CV-selected rank | 19 in 494 of 500 (2, 11, 18 in two each) | 8 |
| pooled relative eigengap, median | 0.330 | 0.744 |
| `angle` null q95, median | 51.0° | 35.4° |
| phase offset in the injected block: median \|d\|, PC1 \|corr\| | 1.126, 0.994 | p180: 1.12, 0.73 |
| phase offset in the other block | 0.031, 0.015 | NMR: 0.21, 0.13 |

The construction matches ADNI on the calibrated quantities (AUROC, the basis following phase, the size of the offset)
and differs on the latent geometry. The CV selects every component it is offered (rank 19 against 8). The stage means
are less one-dimensional (eigengap 0.33 against 0.74), and the `angle` null is wider (51° against 35°). The offset is
the injected block's whole first axis, where ADNI's p180 has other structure of its own.

The paper-grade run has no unadjusted arm. At calibration the unadjusted arm rejected `delta` at 0.15 and `shape` at
0.45 at the anchor, so the construction confounds (criterion 3).

## 2. Size

`size_readings.csv`, bound α + 2·√(α(1 − α)/R) = 0.0695 at R = 500, so at most 34 rejections of 500 hold:

| control | `delta` | `angle` | `shape` |
|---|---|---|---|
| anchor (e = 0) | 17 / 500 = 0.034 | 30 / 500 = 0.060 | 20 / 500 = 0.040 |
| translation (e = 1.0) | 17 / 500 = 0.034 | 27 / 500 = 0.054 | 34 / 500 = 0.068 |

- **Every reading holds.** Translation `shape` holds at the largest count the bound allows (34). Under the rule's
  design (D8), a valid test fails at most one of six readings about 17 % of the time at exact α, so one reading near
  the bound is expected.
- **`delta` is conservative (0.034 at both controls).** Six covariate columns at n = 1,501 cost almost nothing in
  residual df, so the liberal `delta` the SEA-AD configuration showed at n = 80 (`residual_df_shrinkage`) does not
  arise here.
- **The size sentence**, verbatim: "Size holds at the ADNI design: all six readings at or below the one-sided bound".

## 3. Power

`adni_design_rates.csv` (target pairs), SE the binomial standard error, R = 500 at every effect:

| target | e | rate (SE) | CV-selected rank |
|---|---|---|---|
| magnitude / `delta` | 0.1 | 0.034 (0.008) | 19 in 497 |
| | 0.25 | 0.056 (0.010) | 19 in 498 |
| | 0.5 | 0.134 (0.015) | 19 in 500 |
| | 1.0 | 0.530 (0.022) | 19 in 500 |
| orientation / `angle` | 0.25 | 0.070 (0.011) | 19 in 488; 2–3 in 4 |
| | 0.5 | 0.056 (0.010) | 19 in 468; 2–3 in 12 |
| | 0.75 | 0.076 (0.012) | 19 in 418; 2–3 in 53 |
| | 1.0 | 0.058 (0.010) | 19 in 352; 2–3 in 96 |
| shape / `shape` | 0.25 | 0.062 (0.011) | 19 in 487 |
| | 0.5 | 0.048 (0.010) | 19 in 487 |
| | 0.75 | 0.034 (0.008) | 19 in 475 |
| | 1.0 | 0.046 (0.009) | 19 in 472 |

- **Magnitude/`delta`.** It is the only statistic whose rate rises with the effect, but it stays below the floor. At
  e = 1.0 (group B's per-omic stage effect twice group A's) `delta` rejects 0.530. The magnitude surgery also raises
  the stage AUROC, from 0.614 at e = 0.1 to 0.683 at e = 1.0 (`cv_mean_auroc_mean`).
- **Orientation/`angle` and shape/`shape`.** Both sit at the null rate at every effect: 0.056–0.076 and 0.034–0.062.
  These surgeries permute or relocate differential features whose stage effect is δ = 0.02. At an AUROC near 0.60 that
  moves the group's stage means by too little for a test whose null `angle` q95 is about 51°.
- **Rank strata** (`adni_design_rank_strata.csv`; descriptive). The orientation surgery pushes more units to rank 2–3
  as e grows (4 → 96 of 500), but those units do not reject more often. At e = 1.0: 2 / 96 = 0.021 at rank 2–3,
  1 / 8 at 4–7, 26 / 396 = 0.066 at 8+. At e = 0.75: 3 / 53 = 0.057 at 2–3, 35 / 444 = 0.079 at 8+. The Phase 6
  pattern, where power rises above rank 3, does not appear, because there is no power to stratify. Magnitude sits at
  rank 8+ in 499–500 units per effect. The strata's replicates sum to the pooled R in every cell.
- **Off-target rates** (descriptive). They lie at 0.022–0.078: magnitude by `angle` 0.050–0.060 and by `shape`
  0.028–0.042; orientation by `delta` 0.022–0.050 and by `shape` 0.030–0.058; shape by `delta` 0.030–0.044 and by
  `angle` 0.050–0.078. Only one of the 24 off-target readings lies above 0.0695: shape e = 0.75 read by `angle`
  (0.078), which does not enter the verdict.

## 4. Verdict, statements and consequences

**Verdict: `non_rejections_uninformative`.** Size holds at both controls and no target statistic has an informative
effect. The decision file's `statements`, verbatim:

- **`delta`**: "At the ADNI design `delta` reaches at most 0.530 (SE 0.022, R = 500) at e = 1 of the magnitude
  construction, at CV-selected rank 19 (median 19); a non-rejection is uninformative up to the construction maximum
  (e = 1; floor 0.80)."
- **`angle`**: "At the ADNI design `angle` reaches at most 0.076 (SE 0.012, R = 500) at e = 0.75 of the orientation
  construction, at CV-selected rank 2–19 (median 19); a non-rejection is uninformative up to the construction maximum
  (e = 1; floor 0.80)."
- **`shape`**: "At the ADNI design `shape` reaches at most 0.062 (SE 0.011, R = 500) at e = 0.25 of the shape
  construction, at CV-selected rank 2–19 (median 19); a non-rejection is uninformative up to the construction maximum
  (e = 1; floor 0.80)."

**Consequences.**

- **`adni_sexdiff`.** It keeps reading rule 3 in its own change, now as a measurement rather than an absence of one.
  Its `angle` and `shape` non-rejections are uninformative: at this design and stage signal neither statistic detects
  the generator's strongest orientation or shape construction above the null rate. Its `delta` non-rejection (p 0.112)
  is uninformative too, by the same rule: `delta` reaches 0.53 only at a doubled stage effect in one group.
  - Size holds at both controls, and `delta` is conservative there (0.034).
  - So the sensitivity arms' `delta` rejections near α (forced rank 6, p 0.047; `nmr_all_columns`, 0.045) are not
    explained by a liberal null at this design. They remain the arm-setting dependences the case study reads them as.
- **MOTCO.** The documentation states the measured rates beside the SEA-AD statement. Nothing under `src/motco/`
  changes, and the version stays 0.13.0.
- **Roadmap item 10 (Phase 6 exit review).** It reads both case studies as size-controlled designs without power for
  orientation and shape. ADNI also has no useful power for size.

## 5. Limitations and scope

- **Structural analogue.** The blocks are InterSIM omics, not metabolomics: 367 methylation and 131 expression
  features (498) against p180's 133 and NMR's 165 (298). The block shares differ: methylation / expression 74 / 26 %
  against NMR / p180 55 / 45 %. The
  generator's latent geometry differs from ADNI's: the CV rank is at the cap (19 against 8), the eigengap is 0.33
  against 0.74, and the phase offset is the injected block's whole first axis (0.99 against 0.73). No statement is made
  about metabolite-level behaviour.
- **One table.** Only the intersection cohort's table (n = 1,501) is measured. The CN / LMCI / AD arm (n = 1,130) is
  inferred, not measured.
- **One strength.** The run used one calibrated stage signal (δ = 0.02, `s_z` = 0.1, `s_p` = 1.75). The calibration
  matched AUROC, the basis following phase and the offset size. It did not match the rank, which the design records but
  does not tune to.
- **The covariate model.** Covariates enter additively and linearly: one confounded numeric (`z`, weakly injected at
  `s_z` = 0.1), one confounded factor (`category`, declared, not injected), one injected binary nuisance and two noise
  columns.
- **The primary only.** The latent-adjusted, residual-corrected and forced-rank arms, and ADNI's export variants, were
  not run, so their operating characteristics at this design stay inferred.
