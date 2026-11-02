# Phase 6 exit review — 2026-11-02

**Decision: Phase 6 is CLOSED on a size-controlled null; the attribution leg withheld, lifted by a named
demonstration or reported as failed on its own terms.** The exit gate is recorded as **not met as written**. Its
reproducibility criterion is met. Its driver-tracing criterion is **unexercised**: neither case study has a global
trajectory claim to trace, which is different from failed. The phase's purpose is demonstrated up to inference: an
end-to-end, reproducible, size-controlled analysis of two real cohorts whose non-rejections are measured (ADNI) or
bracketed (SEA-AD) to be uninformative. The attribution leg is **withheld**, not waived. Two consequent changes, each
with its own proposal, lift it or report it failed (§6), as the Phase 5 review withheld its magnitude specificity
claim and the [2026-10-08 addendum](phase5-magnitude-remeasurement-2026-10-08.md) lifted it. Phase 7 stays deferred.

**Inputs.** No new measurement was made for this review; it decides on what the following produced.

- MOTCO reports and decision files:
  - the [covariate Type I report](covariate-type-i-2026-10-16.md), `results/covariate-type-i-2026-10-16/covariate_type_i_decision.json` (`covariates_validated`);
  - the [rank-floor report](rank-floor-2026-10-23.md), `results/rank-floor-2026-10-23/delta0.05/rank_floor_decision.json` and `results/rank-floor-2026-10-23/delta0.1/rank_floor_decision.json` (`adopt_floor` at both);
  - the [design-guidance pilot](design-guidance-pilot-2026-10-28.md), its δ arm in `results/design-guidance-pilot-2026-10-28/arms_magnitude/report/` and `results/design-guidance-pilot-2026-10-28/arms_geometry/report/`;
  - the [design-guidance slice](design-guidance-slice-2026-10-29.md), `results/design-guidance-slice-2026-10-29/core_magnitude/report/design_point_operating.csv` and `design_point_stage_signal.csv`;
  - the [case-study Type I report](case-study-type-i-2026-10-29.md), `results/case-study-type-i-2026-10-29/case_study_type_i_overall_decision.json` (`diagnose_configuration_excess`);
  - the [case-study diagnosis](case-study-diagnosis-2026-10-30.md), `results/case-study-diagnosis-2026-10-30/case_study_diagnosis_decision.json` (`residual_df_shrinkage` for the adjusted `delta`, `residualized_basis` for the `stage` arm);
  - the [corrected-permutation report](corrected-permutation-2026-10-30.md), `results/corrected-permutation-2026-10-30/corrected_permutation_decision.json` (`adopt_df`);
  - the [ADNI-design report](adni-design-2026-11-02.md), `results/adni-design-2026-11-02/run/adni_design_decision.json` (`non_rejections_uninformative`, `"indicative": false`);
  - the [Phase 5 exit review](phase5-exit-review-2026-09-10.md), for the form, and the [roadmap](../roadmap.md)'s Phase 6 section.
- The two case-study reports, each read at the commit named (both on `main`, with no tracked-file changes):
  - [Omnia-Lab-UC/sead_sexdiff](https://github.com/Omnia-Lab-UC/sead_sexdiff) at `2d427b05376069f4ad06542f8a9983b00939265c`: `reports/motco-residual-correction-primary-2026-10-31.md` (MOTCO 0.13.0, `--residual-correction df`), cited below as the *SEA-AD report*;
  - [Omnia-Lab-UC/adni_sexdiff](https://github.com/Omnia-Lab-UC/adni_sexdiff) at `1ee24daf2c9ac5a424385cdb622046a53a180a4d`: `reports/motco-adni-2026-11-02-2.md` (MOTCO 0.10.0), cited below as the *ADNI report*.
- The attribution guards on this review's tree: `src/motco/analysis.py:519–523` and `src/motco/stats/attribution.py:368–372`. The first is pinned by `tests/test_analysis.py:579` (`test_attribution_with_covariates_is_rejected_before_fitting`), which passes.

## 1. The exit-gate criteria

Phase 6's purpose is "demonstrate that the validated method answers a meaningful biological question end to end". Its
exit gate has two parts. Each is recorded with one of three values: **met**, **not met** or **unexercised**.

| criterion | status | evidence |
|---|---|---|
| **G1** — collaborators can trace each global trajectory claim to stable feature/module drivers | **unexercised**; attribution **unavailable** at the case-study configuration | No claim exists. SEA-AD rejects no statistic in any of its 44 analyses (SEA-AD report, Result and §11), and ADNI's primary rejects none (ADNI report, Summary). Had a claim existed, attribution could not have run on it: `analyze_trajectories` refuses attribution whenever covariates are declared (`src/motco/analysis.py:519–523`), and both primaries declare six covariate columns (§4, D6). |
| **G2** — collaborators can reproduce the result from documented inputs and configuration | **met** | SEA-AD: the re-issue's gates and reproducibility check are byte-identical on the re-issue machine, and every output that reads no p-value reproduces the superseded outputs made on another machine at rounding level, with every rank and CV AUROC equal (SEA-AD report §12). ADNI: `adni-sexdiff report --check` regenerates the report from the committed outputs and compares it cell by cell (ADNI report, Reproduction). Both run every analysis through one entry point, `motco.analysis.analyze_trajectories`, and record its `motco analyze` CLI equivalent. |
| the gate as a whole | **not met as written** | G2 is met; G1 is unexercised, not met. |
| the purpose | **demonstrated up to inference** | Two real designs are analysed end to end with size measured at both primaries (§2), and their non-rejections are measured or bracketed to be uninformative (§2, §5). The interpretive half, the drivers, has no real-data demonstration (§6). |

**The vacuous pass is refused.** Read literally, "each global trajectory claim" ranges over an empty set, and the gate
would pass with nothing traced. This review does not take that reading. The purpose sentence and the conditional Work
bullet ("For significant orientation, …") both presuppose a detected difference, and a gate that passes because nothing
was found would mean nothing. Nor is G1 recorded as *failed*: failure would say that attribution was tried on a claim
and produced no stable drivers, and nothing of the kind happened.

**The Work list.**

| Work bullet | status | evidence |
|---|---|---|
| Predefine cohorts, stages, covariates, missing-data handling, and preprocessing | **met** | SEA-AD report, reading rules "predeclared in the proposal before any corrected run"; ADNI report, Reading rules (`READING_RULES`, "decided by `adni-analyze` before any run") |
| Fit one pooled, cross-validated PLS representation | **met** | SEA-AD report §3 (astrocyte `cohort_30`: CV rank 5, `cv_mean_auroc` 0.6356); ADNI report, Primary analysis (CV rank 8, `cv_mean_auroc` 0.599514) |
| Report magnitude, orientation, and validated shape results with uncertainty | **met** | RRPP p-values and null summaries in every analysis (SEA-AD report §3; ADNI report, Primary analysis); size measured at both primaries (§2) |
| For significant orientation, produce observed and PLS-reconstructed feature contrasts, bootstrap stability, and pathway/module enrichment | **withheld** — no claim; unavailable with covariates | §4, D1 and D6; §6 |
| Include sensitivity to PLS dimensionality and preprocessing choices | **met** | SEA-AD report §4 (twelve forced-rank analyses) and §6 (sixteen preprocessing-variant analyses); ADNI report, Rank sensitivity (forced ranks 2, 3, 4, 6, 9, 12) and Export sensitivity (six variants) |
| Make the analysis reproducible from raw-input contracts to final figures | **met** | SEA-AD report §1 (pin, gates, commands) and §12; ADNI report, Reproduction |

## 2. The evidence ledger

Every size and power statement carries one of three provenance classes:

- **measured**: a committed decision file at that design and configuration;
- **bracketed**: measured at neighbouring designs whose stage signal or n encloses the cohort's;
- **inferred**: from a mechanism, such as κ = √(n / (n − r)) depending on n and the reduced rank r only.

The paper may state a measured or a bracketed statement, naming its class. It may not turn a bracket into a number for
the cohort. Rates are rejection rates at p ≤ 0.05; `delta` / `angle` / `shape` in that order.

| analysis group | design | result | size | power at its stage signal |
|---|---|---|---|---|
| SEA-AD astrocytes `cohort_30` (primary) | n = 80 (F 11/10/28, M 9/10/12); RNA + ATAC pseudobulk, 5,000 + 5,000 features; six covariate columns; CV rank 5; AUROC 0.6356; `df`, κ 1.06904 | p 0.501 / 0.449 / 0.331; none rejected | **measured** under `df` on this table: pooled 0.031 / 0.0425 / 0.058, every reading at most 32 / 500 (corrected-permutation report §2) | **bracketed**: every rate 0.010–0.062 at n = 80 and anchor AUROC 0.516 / 0.580 (rank-floor report §3–§4); no statistic above 0.09 at n = 200 and AUROC 0.553 / 0.729 (pilot §3.3–§4, δ arm) |
| SEA-AD microglia `cohort_30` (primary cohort of the positive-control cell type) | n = 76; CV rank 2; AUROC 0.6852; κ 1.0731 | p 0.705 / 0.456 / 0.991; none rejected | **measured** under `df` on this table: pooled 0.029 / 0.05 / 0.052 (corrected-permutation report §2) | **bracketed**, as above |
| SEA-AD `cohort_50`, astrocytes and microglia (sensitivity cohorts) | n = 74 (rank 2, AUROC 0.6965, κ 1.0753); n = 63 (rank 4, AUROC 0.6565, κ 1.0903) | p 0.628 / 0.242 / 0.514 and 0.834 / 0.291 / 0.644; none rejected | **inferred** from κ (SEA-AD report §13) | **bracketed**, as above |
| SEA-AD sensitivity arms | 40 analyses: twelve forced ranks, sixteen preprocessing-variant analyses (eight variants × two cohorts), eight latent-adjustment arms, four donor-exclusion analyses | no rejection in any | **inferred** (eight-covariate arms at r = 12, exclusion variants at n = 58–70; SEA-AD report §13); the `latent_adjusted` arm's `angle` and `shape` **inadmissible** (`residualized_basis`) | **bracketed**, as above |
| ADNI primary | n = 1,501 (F 237/310/125, M 223/448/158); p180 + NMR, 133 + 165 features; five covariates in six design columns; CV rank 8; AUROC 0.599514; uncorrected (κ ≈ 1.003); 9,999 permutations | p 0.1116 / 0.6879 / 0.7105; none rejected | **measured**: 0.034–0.068 at both controls against 0.0695 (ADNI-design report §2) | **measured**, `non_rejections_uninformative`: `delta` at most 0.530 (SE 0.022), `angle` 0.056–0.076, `shape` 0.034–0.062 (ADNI-design report §3) |
| ADNI sensitivities | 14 analyses: `latent_stage`; six export variants (`nist_bridge`, `per_run_standardized`, `nmr_all_columns`, `raw_scale`, `log2_plus_one`, `glog`); the CN / LMCI / AD arm (n = 1,130); forced ranks 2, 3, 4, 6, 9, 12 | only `delta` rejects: at forced rank 6 (p 0.0471) and under `nmr_all_columns` (p 0.0447) | **inferred**: the arms were not run at the design (ADNI-design report §5) | measured only for the primary's construction |

The sources of each row's design and result are the SEA-AD report §3, §4, §6, §7, §9 and §13, and the ADNI report's
Summary, Primary analysis and sensitivity sections. The ADNI κ is the ADNI-design report's (§1). The `residualized_basis`
verdict is the [case-study diagnosis](case-study-diagnosis-2026-10-30.md)'s.

**Why SEA-AD's power is a bracket and not a measurement.** The SEA-AD cohorts' held-out AUROCs (0.6356–0.6965) lie
inside the pilot's δ-arm range (0.553–0.729), where no statistic has power at n = 200. They lie above the rank-floor
check's 0.516 / 0.580, where every rate sits at the null at the astrocyte cohort's own table, without covariates. No
run was made at n = 80 and AUROC 0.6356–0.6965 with the case study's covariates. The non-rejections are
uninformative on both sides of the bracket, but the bracket gives no number for the cohort.

**The binding readings**, as the case studies predeclared them:

- **SEA-AD primaries** (SEA-AD report, reading rules and §11). The `angle` and `shape` non-rejections are uninformative
  at this n, and the `angle` null is too wide to inform a non-rejection. The `delta` non-rejections are weak evidence
  against a size difference, since `delta` is conservative under `df` (0.031 / 0.029 pooled). The CV-selected rank,
  eigengap and null summary are reported beside `angle`.
- **SEA-AD arms.** Each is read beside the primary under its own rules. No rank, preprocessing, latent-adjustment or
  exclusion dependence was found. The `latent_adjusted` arm's `angle` and `shape` p-values are reported, not read.
- **ADNI primary** (ADNI report, `READING_RULES`). `delta` is not rejected; the `angle` and `shape` non-rejections are
  uninformative by rule 3, which states that no power statement exists at this design. The ADNI-design check turns
  that absence into a measurement and extends it to `delta`: "a non-rejection is uninformative up to the construction
  maximum" for all three statistics (ADNI-design report §4). `adni_sexdiff` confirms rule 3 in its own change (§8).
- **ADNI sensitivities.** A rejection status that differs from the primary's is dependence on the arm's setting, never
  a replacement of the primary (`READING_RULES` 1), and the rank sweep's rules count no ranks that reject (sweep rules
  2, 5 and 6; §4, D5).

## 3. What Phase 6 establishes, withholds, and does not allow to be claimed

**Claimable now:**

- **An end-to-end, reproducible real-data analysis of two designs through one entry point.** SEA-AD (n = 80, RNA +
  ATAC pseudobulk, 44 analyses) and ADNI (n = 1,501, serum metabolomics, 15 analyses) each ran every analysis through
  `analyze_trajectories`. Each had predeclared reading rules, sensitivity to rank and preprocessing, and a
  reproduction record (§1, G2).
- **Size at both primaries, measured.** SEA-AD `cohort_30` under `df`: pooled 0.029–0.058 per cohort table and
  statistic against 0.0597 (corrected-permutation report §2). ADNI: 0.034–0.068 against 0.0695 (ADNI-design report §2).
- **The null result, stated with its uninformativeness.** Neither primary rejects any statistic: SEA-AD astrocyte and
  microglia `cohort_30`, and ADNI. The non-rejections are uninformative. At ADNI this is measured: no statistic reaches
  the 0.80 floor at any construction effect. At SEA-AD it is bracketed, and stated as a bracket (§2).
- **The design finding** (§5): at these designs the held-out stage signal, not n, limits power.
- **The configuration lessons.** Each is a measured MOTCO result that a real analysis ran into first:
  - where the stage signal is weak, the CV rank needs a floor of 2 (rank-floor report, `adopt_floor`; MOTCO 0.10.0);
  - with six covariate columns at n = 80, the Freedman–Lane null shrinks with the reduced model's residual df, and the
    `df` correction restores size (case-study diagnosis, `residual_df_shrinkage`; corrected-permutation report,
    `adopt_df`; MOTCO 0.13.0);
  - the covariate-blind basis follows technical axes where they dominate the variance (SEA-AD library composition,
    ADNI assay phase). The axis diagnostics report them, and no statistic is re-read for them. The `stage` latent
    adjustment is inadmissible for `angle` and `shape` at n = 80 (`residualized_basis`);
  - the held-out CV stage AUROC is the observable that places a cohort on the power axis (§5).

**Withheld, with a named path to lift it:**

- **Driver tracing on real data**: G1, the Phase 6 Work bullet for attribution, and the roadmap's objective 4. It is
  withheld because it was never exercised. It has not failed, and it is not waived. Items 11 and 12 lift it, or report
  it failed on its own terms (§6). The paper's real-data section states it as withheld and names the gap.

**Not claimable:**

- **Any agreement reading.** A non-rejection of `delta`, `angle` or `shape` in either cohort is not evidence that the
  sexes' trajectories agree in size, direction or shape. Both case studies' rules say so, and the measured or
  bracketed power says why.
- **The sensitivity-arm rejections as results.** ADNI's `delta` rejections at forced rank 6 (p 0.0471) and under
  `nmr_all_columns` (p 0.0447) are setting dependence under the case study's predeclared rules (§4, D5).
- **A numeric power for SEA-AD.** The bracket says that the non-rejections are uninformative. It gives no power
  number at n = 80 and AUROC 0.6356 (§2).
- **Metabolite- or feature-level claims from the structural analogue.** Both operating-characteristic constructions
  use InterSIM blocks (498 features), and the ADNI one differs from ADNI on its latent geometry (§4, D7). They support
  statements about what the analysis can detect at the design. They support nothing about metabolite-, gene- or
  peak-level behaviour.

## 4. Resolution of each deviation

**D1 — No global trajectory claim to trace. Claim revision.**
Neither primary rejects: SEA-AD astrocyte `cohort_30` reads p 0.501 / 0.449 / 0.331 and ADNI p 0.1116 / 0.6879 /
0.7105. SEA-AD rejects nothing in any of its 44 analyses. The mechanism is power, not size: size is measured at both
primaries, and power is absent at both designs (§2, §5). G1 is recorded as unexercised and the phase closes on the
null. No further analysis of either cohort is run to produce a claim. The case studies' predeclared rules forbid that.
SEA-AD's reading rules state that "nothing is re-tuned", and ADNI's `READING_RULES` 9 reads "The analysis is run once at
these settings; a second run at other settings is a new proposal". Sources: SEA-AD report, reading rules and §11;
ADNI report, Summary and Reading rules.

**D2 — The rank-1 `angle` degeneracy. Method revision (MOTCO 0.10.0).**
Under the 0.9.0 rule both SEA-AD cohorts selected rank 1. There `angle` is 0° or 180° by construction, and its RRPP
null is a two-point distribution. `minimum-latent-rank` floors the CV candidates at 2. The rank-floor check returned
`adopt_floor` at δ = 0.05 and 0.1: the floor-2 arm holds 0.0695 at the anchor (`delta` 0.042 / 0.044, `angle` 0.054 /
0.060, `shape` 0.036 / 0.038). The floor defines `angle` and `shape` in at least two dimensions; it does not make them
powerful. Source: [rank-floor report](rank-floor-2026-10-23.md) §1–§4;
`results/rank-floor-2026-10-23/delta0.05/rank_floor_decision.json`,
`results/rank-floor-2026-10-23/delta0.1/rank_floor_decision.json`.

**D3 — Residual-df shrinkage of the Freedman–Lane null. Method revision (MOTCO 0.13.0); the case study was re-issued
with no status change.**
At the case study's configuration (n = 80, six covariate columns) the covariate-adjusted primary exceeded the bound on
translation `delta` at δ = 0.1: 0.072 against 0.0695 ([case-study Type I report](case-study-type-i-2026-10-29.md),
`diagnose_configuration_excess`). The [diagnosis](case-study-diagnosis-2026-10-30.md) found `residual_df_shrinkage`.
The reduced model's residuals have mean square RSS / n rather than RSS / (n − p_red), with p_red = 10. Pooled `delta`
rises 0.0545 → 0.0595 → 0.066 with one, three and six covariate columns. The opt-in `df` correction multiplies the
residuals by κ = √(n / (n − r)). The [corrected-permutation check](corrected-permutation-2026-10-30.md) returned
`adopt_df`, with pooled rates of 0.029–0.058. SEA-AD then re-issued every MOTCO-dependent output under `df`: every
p-value moved up or stayed, and no rejection status changed in any of the 44 analyses (SEA-AD report, Result and §12).
ADNI keeps the uncorrected null, because κ ≈ 1.003 at n = 1,501 is immaterial (ADNI-design report §1) and `delta` is
conservative there (0.034, §2).

**D4 — The covariate-blind basis follows technical axes, and the `stage` latent adjustment is inadmissible for `angle`
and `shape`. Claim revision.**
SEA-AD: the first latent component of the astrocyte `cohort_30` primary tracks `ribosomal_mean` (|ρ| 0.840858) more
closely than stage (ρ 0.491212). In microglia it tracks `mitochondrial_pseudobulk` (0.843828, against stage 0.688175)
(SEA-AD report §5). The positive control, microglia against astrocytes, is a control for the stage signal. It does not
hold on the primary basis in `cohort_30`, through E3 alone (0.6882 against `mitochondrial_pseudobulk` 0.8438), and it
holds on the `latent_adjusted` basis (E3 0.5279 against `tmm_rna` 0.4153) (SEA-AD report §10). ADNI: the primary's
first component tracks the assay phase (|ρ| 0.542; ADNI report, Axis diagnostics of the primary). The `stage` latent
adjustment that would remove such axes is `residualized_basis` at n = 80: at a forced rank of 2 it still exceeds the
bound on anchor `angle` (0.084 / 0.082; case-study diagnosis). The claim is revised, not the method:
- an axis correlated with a technical variable is a measurement-space finding, reported beside the statistics and
  never a reinterpretation of them;
- the `latent_adjusted` arm's `angle` and `shape` p-values are reported, not read;
- the covariate-blind analysis stays the primary in both case studies.

**D5 — ADNI's two sensitivity-arm `delta` rejections. Claim revision: none; the predeclared reading stands.**
`delta` rejects at forced rank 6 (p 0.0471), one of six forced ranks. It also rejects under the `nmr_all_columns`
export variant (p 0.0447), one of six export variants (ADNI report, Rank sensitivity and Export sensitivity). Under
`adni_sexdiff`'s predeclared sweep and membership rules both are setting dependence. The primary's non-rejection
stays the result, and the sweep's rules count no ranks that reject (sweep rules 2, 5 and 6).
- *A liberal null does not explain them.* Size holds at the design, and `delta` is conservative there (0.034;
  ADNI-design report §2).
- *They are not a claim.* No multiplicity-controlled rule makes one.
- *They are not a G1 candidate.* Running attribution on them would be the forking path this review closes, and it
  would meet the gap of D6.

**D6 — Attribution is unavailable with covariates. Capability gap → item 11.**
Two guards block it:
- `analyze_trajectories` raises "Attribution is not available with covariates: it works from unadjusted cell means"
  whenever covariates are declared (`src/motco/analysis.py:519–523`, pinned by `tests/test_analysis.py:579`).
- The lower-level `analyze_orientation_attribution` accepts covariate-adjusted means only as a precomputed
  `mean_table`, and then requires `bootstrap_replicates = 0`, because "row-level bootstrap semantics are unavailable"
  (`src/motco/stats/attribution.py:368–372`).

Both case-study primaries declare six covariate columns. So a claim from either primary could not have been traced to
*stable* drivers as the analysis is built today. Unadjusted attribution is no substitute: under a confounded covariate
the unadjusted contrast is a different estimand, and the covariate Type I check's unadjusted arm rejected `shape` at
0.364 on null data ([covariate Type I report](covariate-type-i-2026-10-16.md) §2).

**D7 — The structural analogue is wider at ADNI. Claim scope.**
The ADNI-design construction matches ADNI on the calibrated quantities: held-out AUROC 0.606 against 0.600, and
first-component |ρ| with phase 0.578 against 0.542. It does not match ADNI's latent geometry:
- CV rank 19 in 494 of 500 anchor units, against 8;
- eigengap 0.330 against 0.744;
- `angle` null q95 51.0° against 35.4°.

Its blocks are InterSIM omics (498 features), not ADNI's 298 metabolites (ADNI-design report §1, §5). Its measured
rates say what the analysis can detect at the design; they say nothing about metabolite-level behaviour. The SEA-AD
statements carry the same limitation (the same 498 InterSIM features). The ADNI-design report records it as wider at
ADNI, because there the construction's rank sits at the cap and its eigengap is less than half of ADNI's.

## 5. The cross-study finding and the binding reading

**At these designs the held-out stage signal, not n, limits power.**

| run | n | held-out stage AUROC | best power | source |
|---|---|---|---|---|
| ADNI-design check | 1,501 | 0.606 (ADNI 0.600) | `delta` at most 0.530 (SE 0.022), at a doubled stage effect in one group (e = 1.0); `angle` and `shape` at the null rate | [ADNI-design report](adni-design-2026-11-02.md) §1, §3 |
| design-guidance slice (three stages, two blocks) | 80 | 1.000 (mean ≥ 0.9998) | `delta` 0.620 (SE 0.022) at e = 0.1 | [slice report](design-guidance-slice-2026-10-29.md) §3.1, §4 |
| design-guidance pilot, δ arm (δ = 0.05 / 0.1) | 200 | 0.553 / 0.729 | no statistic above 0.09 at any effect (magnitude ≤ 0.03, orientation ≤ 0.07, shape ≤ 0.09; R = 100) | [pilot report](design-guidance-pilot-2026-10-28.md) §3.3, §3.4, §4 |

A cohort almost nineteen times larger, at AUROC 0.60, detects a doubled stage effect less often than n = 80 at a
saturated stage signal detects a 10 % size change. The SEA-AD cohorts (AUROC 0.6356–0.6965) sit inside the δ arm's
bracket, at a smaller n than it measured.

**Binding reading**, for the paper and for any later real-data use of `motco analyze`: *a non-rejection is informative
only where a measured operating characteristic at the cohort's own design and stage signal says so; the held-out CV
stage AUROC is reported beside the CV-selected rank with every result*.
- It extends the Phase 6 statement that already requires the rank beside any `angle` result, and replaces nothing in
  it.
- It is a reading condition, not a power model. Three regimes do not fit a power curve in AUROC, and none is offered.

## 6. The withheld leg and its lifting criterion

**The gap.** The lifting demonstration must run at the configuration in which a real claim would be read: covariates
declared, the rank floor, net displacement and, at small n, `df`. Attribution as built cannot run there (D6: the
guards at `src/motco/analysis.py:519–523` and `src/motco/stats/attribution.py:368–372`). Two consequent changes follow,
in order, each with its own proposal.

**Item 11 — covariate-adjusted attribution in `motco analyze`** (MOTCO). It has three parts:
- attribution from covariate-adjusted cell means, evaluated at the record's reference row;
- bootstrap stability that refits the adjustment within each resample on the frozen basis;
- its own validation that the attribution of an adjusted contrast is stable and unbiased under a confounded covariate.

It is byte-identical without covariates. It would modify `real-data-trajectory-analysis` and
`orientation-driver-attribution`. Its proposal decides, after a smoke timing, whether the validation needs cluster
compute.

**Item 12 — the plasmode demonstration** (after item 11). It lifts the leg or reports it failed.
- *Data.* The real exported blocks, covariates and design of SEA-AD astrocyte `cohort_30`, the primary cohort. SEA-AD
  is open data, while ADNI's are under a data-use agreement. This configuration's size under `df` is measured. And
  ADNI's analogue sits at the CV cap, so its latent geometry is the less understood of the two.
- *Construction.* A planted group × stage difference of a declared kind, with a known feature set, is added to the
  real data at declared strengths. The orientation kind comes first, because the Work bullet is conditioned on
  orientation. Each plasmode unit is analysed by `motco analyze` at the case-study configuration with attribution.
- *Lift criteria.* The leg is lifted when all three of these hold:
  - (i) the unplanted construction holds size at the one-sided bound α + 2·√(α(1 − α)/R);
  - (ii) the planted difference is detected at a rate fixed in item 12's proposal;
  - (iii) the planted features are recovered among the attribution's top drivers, with bootstrap sign stability and
    top-k frequency above thresholds fixed in that proposal before any run.

  Otherwise the leg is **reported failed on its own terms**, with the failing criterion named.
- *What this review leaves to item 12's proposal.* Three things are left open. The first is the detection rate in
  (ii). The second is the thresholds in (iii), which are fixed there before any run. The third is the null-making
  step for (i): real sex labels may carry a real difference, so "unplanted" data need one, such as a within-stage
  label permutation that preserves the covariates. The proposal may also add an unadjusted arm as a descriptive
  comparison. That arm cannot lift the leg (D6).
- *Ownership.* MOTCO owns the planting construction and the decision rule, as it does for its other validation
  checks. It reads `sead_sexdiff`'s export through a documented input contract and commits aggregates only.
  `sead_sexdiff` supplies the export, and no donor-level data are committed to MOTCO. Item 12's proposal confirms the
  split.

If item 11 proves harder than expected, the leg stays withheld and the paper reports it so, with the gap named. The
closure does not depend on it.

## 7. Phase 7 stays deferred — unchanged

SNF-native trajectory statistics stay deferred on the [Phase 5 exit review](phase5-exit-review-2026-09-10.md)'s §5
grounds. Nothing in either case study bears on SNF: both measured the PLS latent space throughout.

## 8. Consequent work, in order

1. **Item 11 — covariate-adjusted attribution** (MOTCO; owed by the withheld leg; §6).
2. **Item 12 — the plasmode demonstration** (MOTCO, with `sead_sexdiff`'s export; owed by the withheld leg; §6). It
   starts after item 11.
3. **`adni_sexdiff` reading rule 3** (owed by that repository, in its own change). It confirms the rule, now citing
   `results/adni-design-2026-11-02/run/adni_design_decision.json`, so that its `angle`, `shape` and `delta`
   non-rejections are read as measured to be uninformative at its design.
4. **Independent of Phase 6; none of these blocks anything Phase 6 owes:**
   - item 5, optional block weighting;
   - item 6, defaulting the null-dispersion materiality rule;
   - item 7(b), group-specific dispersion in the generator and its null study;
   - item 8's follow-ups: whether `df` becomes the analysis default, and the `cohort_50` tables, whose size under `df`
     is inferred and on which no rejection depends;
   - the toy-fixture decision (cross-cutting engineering work).

`sead_sexdiff` and `adni_sexdiff` cite this review in their reports' limitations when they next re-issue. Nothing in
either changes now.

## 9. What this review does not do

- It makes no new measurement and changes no code. Nothing under `src/motco/` or `tests/` moves, and no release is
  made.
- It edits no predeclared artefact: no case-study reading rule, decision file or committed report. It does not
  rewrite the Phase 6 gate's text, which stays as written beside the status recorded here.
- It re-reads no case-study result under new rules. Every reading above is the one the case studies predeclared,
  confirmed or cited.
- It ends the two cohorts' use as Phase 6 evidence. No further analysis of SEA-AD or ADNI is read as Phase 6
  evidence, so a later analysis of either is a new study with its own proposal, not a continuation of this phase.
  The plasmode of item 12 uses SEA-AD's real blocks as the substrate for a planted signal. It reads no sex × stage
  result from them.
