# Covariate-adjusted attribution at the case study's configuration — 2026-11-02

> This report validates covariate-adjusted orientation attribution as `motco analyze` runs it with declared
> covariates (openspec change `covariate-adjusted-attribution`, spec `adjusted-attribution-validation`; roadmap Next
> changes item 11). Phase 6 closed with its attribution leg withheld because attribution could not run on a
> covariate-adjusted contrast ([exit review](phase6-exit-review-2026-11-02.md) §4 D6). MOTCO 0.14.0 adds it. This check
> asks whether it removes a confounded covariate from the drivers, and whether it recovers and stabilizes planted
> drivers as well as an oracle that knows the injection.
>
> Every number below is read from `results/adjusted-attribution-2026-11-02/`:
> - `run/adjusted_attribution_decision.json`, `run/adjusted_attribution_summary.csv`, `run/merge.json`,
>   `run/decide.json`, `run/settings.json`, `run/sacct.txt`, `run/PROVENANCE.txt`;
> - `calibration/calibration.json` (every calibration iteration) and `calibration/probe/` (the diagnostic behind the
>   gate amendment);
> - the smoke `smoke/` (R = 20 × 99 permutations × B = 50, indicative) and `timing/timing.json`.
>
> `run/records/`, `run/merged.jsonl` and the 5.4 MB `run/adjusted_attribution_metrics.csv` stay on `ing`;
> `run/PROVENANCE.txt` names them with their sha256. Readings are means over R = 500 replicates per cell at 999
> permutations and B = 200 bootstrap replicates, on the **observed** component of the **principal-orientation** block
> unless stated otherwise.

**Run.**
- One SLURM array on `ing` (job 932521, partition `512x1024`, n2 and n11 excluded, 50 shards of 30 units), from the
  clean worktree of `feat/covariate-adjusted-attribution` at `3862c6c`. Its git tree `feaf93b8…` equals the
  workstation's, and the change's tests passed there (78 passed).
- All 50 tasks ended `COMPLETED 0:0` between 20:40Z and 20:53Z on 2026-11-02, after a one-unit timing task (job
  932517). 1,500 of 1,500 units, zero consistency failures, zero rank-deficient bootstrap replicates.
- 10.44 elapsed core-hours (planned about 3.5 serial-equivalent; all 50 tasks shared one node); the units' own timers
  sum to 8.42 h.
- `merge` and `decide` ran on `ing`. The decision file's settings hash `5a84844d…` equals `settings.json`'s.

**Verdict: `adjusted_attribution_validated`**, read verbatim from `adjusted_attribution_decision.json`
(`"indicative": false`, `"failing": []`, 500 replicates per cell, adopted calibration iteration 9).

**In one paragraph.**
- **The construction.** The case study's configuration: the astrocyte `cohort_30` table (n = 80), methylation +
  expression, a covariate `z` confounded with group × stage and injected into 20 % of the features at `s_z` = 10, five
  noise covariates (six design columns), rank floor 2, net displacement and the `df` correction. The stage effect is
  the Phase 6 one (δ = 2.0), so there is an orientation to plant. Each unit ran the shipped analysis once, then
  attributed its frozen model three ways: adjusted (the analysis's own attribution), an oracle with the injection
  removed exactly, and unadjusted.
- **The confounding is real.** Without adjustment the covariate takes 0.43 of the top-20 at the anchor beyond the
  oracle, and tilts the contrast toward the injection direction by about 0.51 in every cell.
- **Adjustment removes it.** The adjusted arm's contamination is within 0.005 of the oracle's, and its residual
  projection is −0.003 to −0.002, under 1 % of the unadjusted arm's 0.51–0.52.
- **At no cost to recovery or stability.** Against the oracle, truth precision is +0.038 / +0.005, sign stability
  +0.034 / +0.017 and top-k frequency −0.011 / −0.009 at e = 0.5 / 1.0, all inside their margins.
- **The gates were amended before any unit.** Calibration under the declared rule could not be satisfied: at the
  orientation cells the planted drivers crowd every injected feature out of the unadjusted top-20 at any injection
  strength. With the user's approval the contamination gate moved to the anchor and a projection gate was added at
  every cell (§1); every criterion and margin is as declared.
- **Consequence.** Item 12, the plasmode demonstration on SEA-AD's real astrocyte blocks, may start on 0.14.0.

## 1. Construction

`scripts/adjusted_attribution_check.py` composes the case-study construction by import (`case_study_type_i_check`,
built on the D10 construction of `covariate_type_i_check`):

- **Generator.** The semi-synthetic generator at the astrocyte table `[[11, 10, 28], [9, 10, 12]]`, three stages,
  methylation (367 CpGs) and expression (131 genes), joint magnitude, `p_dmp` 0.1 and the Phase 6 per-omic stage
  effect δ = 2.0. The cells are the zero-effect anchor and the `orientation` surgery at e = 0.5 and 1.0.
- **Covariates, in model order** (six numeric design columns): `z = a_g·g + a_s·s + a_gs·g·s + ε` (D10's
  coefficients 0.5 / 1.0 / 2.0, noise SD 1), injected as `X_f += s_z·sd_f·w_f·z` on a fixed 20 % feature subset
  (100 of 498 features); and `x1` … `x5`, i.i.d. N(0, 1), never injected.
- **Seeds.** Every seed derives from base seed 2033 (no committed check uses it) and the replicate index, so the three
  cells of a replicate share the baseline draw, the covariates, the injection, the RRPP seed and the bootstrap seed.
- **Truth.** The generator's driver set: features whose group-stage differential mean change differs between groups
  (the union over transitions for the principal block). It holds 156 features on average at e = 0.5 and 291 at
  e = 1.0, none at the anchor; 32 and 59 of the 100 injected features are drivers there.

**Arms, on one measurement per unit.**
- **`adjusted`.** `analyze_trajectories` with the six covariates, attribution (B = 200, top-k 20), `min_components`
  2, net displacement, `residual_correction="df"`, the harness CV knobs and 999 permutations.
- **`oracle`.** The attribution primitive on `X̃ − z ⊗ (v / σ)`: the analysed features minus the exact raw-unit
  injection `v = s_z·sd·w` scaled by the pooled analysed scale `σ`, so its cell means are the uninjected ones up to a
  constant vector. Same frozen estimator, bootstrap seed, B and top-k.
- **`unadjusted`.** The primitive on the analysed features as they are (descriptive).
- All three draw the same within-cell bootstrap resamples. Each unit first checks that the exposed estimator
  reproduces the latent scores (to 1e-10 of the largest) and that the adjusted means project to the analysis's stage
  means (to 1e-8 relative); no unit failed either check.

**The measured analyses.** The held-out CV stage AUROC is 1.000 in every cell. The CV-selected rank is 2–3 at the
anchor (292 / 208 units), 2–4 at e = 0.5 (median 3) and 2–6 at e = 1.0 (median 3). The adjusted analysis's `angle`
rejects at 0.046 at the anchor, 0.256 at e = 0.5 and 0.666 at e = 1.0.

### Calibration and the gate amendment

Calibration ran R = 20 × 99 permutations × B = 50 at the anchor and orientation e = 1.0, starting at `s_z` = 10 and
stepping by 5. Every iteration is in `calibration/calibration.json`.

**Under the declared rule** (contamination gap, unadjusted − oracle, at least 0.25 in both cells), iterations 0–8 at
`s_z` 10–50 ended `stop_exhausted`:

| `s_z` | contamination gap, anchor | contamination gap, e = 1.0 | oracle truth precision, e = 1.0 |
|---|---|---|---|
| 10 | 0.485 | 0.000 | 0.975 |
| 15 | 0.4225 | 0.000 | 0.975 |
| 20 | 0.3525 | 0.000 | 0.970 |
| 25 | 0.280 | 0.000 | 0.970 |
| 30 | 0.2675 | 0.000 | 0.9725 |
| 35 | 0.2375 | 0.000 | 0.975 |
| 40 | 0.205 | 0.000 | 0.9775 |
| 45 | 0.1825 | 0.000 | 0.985 |
| 50 | 0.165 | 0.000 | 0.985 |

Raising `s_z` lowered the anchor's gap and never opened one at e = 1.0: in standardized units an injected feature
converges to `z`'s own cell pattern as `s_z` grows. The diagnostic in `calibration/probe/` (six replicates per row,
no permutations, indicative) showed why. At e = 1.0 only about 39 injected features are not drivers, and their best
rank in the unadjusted contrast was 30th to 90th. Weakening the stage effect to δ = 0.25 did not change that (0.0–0.07).
The unadjusted − oracle projection on the injection direction was about 0.5 in every cell nonetheless. The confounding
reaches the contrast; it does not reach the top-20 when a strong orientation is planted.

**The amendment** (2026-11-02, approved by the user before any smoke or paper-grade unit; `design.md` D10). The
contamination gate reads the anchor only, and a projection gate, |unadjusted − oracle projection| − 2·SE ≥ 0.2, reads
every cell. Every criterion and margin is unchanged. Under this rule (version 2), iteration 9 at `s_z` = 10 met every
criterion and was adopted: contamination gap 0.485 at the anchor, projection gaps 0.563 and 0.547, oracle precision
0.975. Its readings equal iteration 0's.

## 2. Gates

| gate | cell | mean | SE | reading | bound | verdict |
|---|---|---|---|---|---|---|
| contamination, unadjusted − oracle − 2·SE | anchor | 0.4255 | 0.0115 | 0.4025 | ≥ 0.25 | holds |
| \|projection, unadjusted − oracle\| − 2·SE | anchor | 0.5091 | 0.0057 | 0.4978 | ≥ 0.2 | holds |
| | e = 0.5 | 0.5207 | 0.0057 | 0.5093 | ≥ 0.2 | holds |
| | e = 1.0 | 0.5113 | 0.0056 | 0.5002 | ≥ 0.2 | holds |
| oracle truth precision | e = 1.0 | 0.9845 | — | 0.9845 | ≥ 0.5 | holds |

The construction is informative. The covariate fills the unadjusted top-20 where nothing else competes (0.475
against the oracle's 0.050 at the anchor), and it bends the unadjusted contrast in every cell. The planted drivers are
recoverable: the oracle's top-20 at e = 1.0 is 98.5 % drivers.

## 3. Criteria

All readings are adjusted − oracle, paired over the cell's 500 replicates, with 2 paired SEs.

| criterion | cell | mean | SE | reading | bound | verdict |
|---|---|---|---|---|---|---|
| contamination (+2·SE) | anchor | 0.0008 | 0.0028 | 0.0064 | ≤ 0.05 | holds |
| | e = 0.5 | −0.0050 | 0.0013 | −0.0024 | ≤ 0.05 | holds |
| | e = 1.0 | −0.0003 | 0.0003 | 0.0003 | ≤ 0.05 | holds |
| residual confounding (\|mean\| + 2·SE) | anchor | −0.0032 | 0.0012 | 0.0056 | ≤ 0.0509 | holds |
| | e = 0.5 | −0.0025 | 0.0011 | 0.0047 | ≤ 0.0521 | holds |
| | e = 1.0 | −0.0022 | 0.0011 | 0.0044 | ≤ 0.0511 | holds |
| truth precision (−2·SE) | e = 0.5 | 0.0384 | 0.0135 | 0.0115 | ≥ −0.05 | holds |
| | e = 1.0 | 0.0047 | 0.0019 | 0.0009 | ≥ −0.05 | holds |
| sign stability (−2·SE) | e = 0.5 | 0.0341 | 0.0039 | 0.0263 | ≥ −0.05 | holds |
| | e = 1.0 | 0.0166 | 0.0040 | 0.0087 | ≥ −0.05 | holds |
| top-k frequency (−2·SE) | e = 0.5 | −0.0110 | 0.0050 | −0.0210 | ≥ −0.10 | holds |
| | e = 1.0 | −0.0087 | 0.0044 | −0.0174 | ≥ −0.10 | holds |

The residual-confounding bound is 0.10 × the unadjusted arm's |mean projection difference| in the same cell (0.5091,
0.5207, 0.5113), so adjustment removes over 99 % of the confounding the unadjusted contrast carries. The adjusted
arm's truth precision and sign stability sit slightly above the oracle's, plausibly because the adjustment also
removes the features' chance within-cell association with the five noise covariates, which the oracle keeps; the check
does not test that explanation. Top-k frequency is 0.009–0.011
below the oracle's, a tenth of its margin: the six refitted coefficients per replicate add some spread to the
per-feature ranks.

The adjusted contrast is not the oracle's contrast. Their cosine is 0.55 at the anchor and 0.74 / 0.79 at e = 0.5 /
1.0, and their top-20 overlap 0.28 / 0.50 / 0.50; the unadjusted arm's are lower (0.49 / 0.69 / 0.74 and 0.17 / 0.43
/ 0.44). Estimating the coefficients on n = 80 moves the per-feature contrast; the criteria ask whether that movement
reintroduces the confounder or loses drivers, and it does neither.

## 4. Verdict and consequence

**Verdict: `adjusted_attribution_validated`.** Both confounding gates and the recoverability gate hold, and every
criterion holds at every cell. The decision file's consequence, verbatim: "covariate-adjusted attribution is validated
at the case study's configuration; item 12 (the plasmode demonstration on SEA-AD astrocyte cohort_30) may start on the
released version (0.14.0)".

- **Item 12** is proposed next, pinned to MOTCO 0.14.0. Its proposal fixes the detection rate, the recovery thresholds
  and the null-making step on the real blocks. The Phase 6 attribution leg stays withheld until it reports.
- **`motco analyze`.** With declared covariates, `--attribution` reports covariate-adjusted drivers; the record's
  `attribution.adjustment` block names the means it used. The validation status in `docs/api/analysis.md` now cites
  this report.

## 5. Descriptive readings

None of these enters the verdict.

- **The unadjusted arm.** Read against the same bounds, it fails contamination at the anchor (0.4485 against 0.05) and
  residual confounding everywhere (0.52–0.53 against about 0.05). Its truth precision is *higher* than the oracle's
  (+0.084 / +0.013): injected features that are also drivers get a boost and are counted as correct, which is why the
  contamination metric counts only injected non-drivers.
- **The PLS-captured component.** Contamination holds (≤ 0.019). Residual confounding reads 0.045 / 0.037 / 0.023
  against 0.019 / 0.019 / 0.031. On the means, adjustment removes 85 / 89 / 97 % of the captured component's
  confounding, short of 90 % at the anchor and e = 0.5.
  The captured view projects the adjusted means through a basis fitted without the covariates, so a share of the
  injection direction that lies in that basis is reconstructed with the signal. The driver report uses the observed
  component (design D9).
- **The residual component.** Contamination holds. On the means, adjustment removes 90 / 93 / 94 % of its
  confounding, but the 2-SE readings, 0.049 / 0.040 / 0.025, exceed the bounds 0.032 / 0.033 / 0.020; truth precision
  at e = 0.5 reads −0.052 against −0.05. The residual is observed minus captured, so it inherits the captured view's
  remainder with the opposite sign.
- **The transitions.** Both transition blocks hold every criterion at every cell, with margins wider than the
  principal block's (residual confounding 0.003 against about 0.028).
- **The `angle`-rejected units.** Restricted to units whose adjusted analysis rejects `angle` at 0.05 (23 / 128 / 333
  units), every criterion holds (largest contamination reading 0.0137; truth precision at e = 0.5 −0.0131).
- **Rank-deficient replicates.** None in any unit. All six covariates are numeric, so no resample can drop a factor
  level.

## 6. Limitations and scope

- **Structural analogue.** The blocks are InterSIM omics (367 methylation and 131 expression features, 498), not
  SEA-AD's 10,000 RNA and ATAC features. Item 12 runs on the real blocks.
- **One table, one stage effect.** Only the astrocyte table (item 12's cohort) at δ = 2.0 is measured. At SEA-AD's own
  stage signal (held-out AUROC 0.63–0.70) the stage transitions are near zero and every arm's contrast is noise; how
  attribution behaves there is item 12's question.
- **Numeric covariates only.** The rank-deficient path for a rare factor level is unit-tested, not measured here.
- **The covariate model.** Covariates enter additively and linearly, as in the trajectory models; the interpretation
  boundary in every adjusted result says so. Adjustment removes only declared covariates.
- **The amended gates.** The contamination gate is read at the anchor only. At the orientation cells, removal is judged
  by the projection, which does not depend on the top-k, while contamination there is a criterion without a gate.
- **The observed principal block.** The verdict reads one component and one block. The PLS-captured and residual views
  hold contamination but not the 90 % residual-confounding bound at every cell (§5).
