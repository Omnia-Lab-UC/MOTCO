# The ADNI latent plasmode — report (2026-11-18)

**Change:** `adni-latent-plasmode-capture` (MOTCO 0.25.0, no package change). **Reads:**
`results/adni-latent-plasmode-2026-11-18/` (`adni_latent_plasmode_decision.json`; per profile `settings.json`,
`preflight.json`, `capture_table.csv`, `point_table.csv`, `accounting.json`; `cohort/reproduction.json`;
`calibration/calibration_table.csv`; `PROVENANCE.txt`). **Parents:** the
[ADNI real-noise anchor](adni-real-noise-anchor-2026-11-17.md) (its point table, bound by sha256) and
[ADNI gain and specificity](adni-gain-specificity-2026-11-17.md) (its plants cell table). **Concept page:**
[measurement visibility](../measurement-visibility.md#the-latent-plasmode-on-adni-gate-failed).

How much of a *real* stage program production's fitted span sees is not established (roadmap open item 2). The only
capture of a weak program is simulated, 0.10 at SEA-AD's stage signal (item 14). Item 22 read the share of a held-out
configuration inside the span fitted on 750 ADNI participants, 0.536, but that configuration is mostly estimation
noise. This change estimated the capture of the program itself from two independent held-out halves: with `m1`, `m2`
the halves' centred covariate-adjusted stage configurations and `Q` a space's basis fitted on a third, disjoint
sample, `tr(m1·Q·Qᵀ·m2ᵀ)` and `tr(m1·m2ᵀ)` are unbiased for the program's energy inside the space and in total. A
calibration on synthetic cohorts with the truth known gated the estimator before any reading.

**In one paragraph.** **The verdict is `gate_failed`, so nothing is read.**

- **The run is clean.** 2,000 cohort units and 7,200 calibration units completed with none failed or
  inconsistent. Both pre-flights passed, and both reproductions hold exactly.
- **The gate fails on bias alone.** `pls_cv`'s total capture misses the ±0.05 limit at five of eight
  strength × point cells. At strength 1 the bias is small and negative (−0.035 to −0.077). At 0.67 it reads −0.45 to
  +0.11, but there the 24-cohort mean has a standard error of 0.12 to 0.37, and an exactly unbiased estimator would
  have passed all four points with probability about 0.001.
- **It would not have resolved either.** The program is not present at 80 participants under the rule (energy
  0.779 against a threshold of 1.015), and no point's interval is narrower than ±1.13 against the limit ±0.15.
- **The consequence fixed with the rule.** A diagnosis of the estimator comes before any reading. No real-noise
  capture is stated, the simulated 0.10 stays the only capture of a weak program, and item 22's held-out readings stay
  mixtures of program and noise.

## 1. Configuration and provenance

- **Export.** `adni_sexdiff/results/motco_input` on ing, read in place through `ADNI_CONTRACT` and never copied: 1,501
  participants, the `p180` (133) and `nmr` (165) blocks, Female 237 / 310 / 125 and Male 223 / 448 / 158 at
  CN / MCI / AD. Its five file hashes are in both settings and were unchanged after the run.
- **Profile `cohort`.** The anchor driver's units, unchanged: 500 per point at 80, 200, 500 and 750 participants
  subsampled within sex × stage cells. The units at 80, 200 and 500 are the anchor run's, those at 750 the gain run's.
  Each unit fits the spaces on its sample and splits the remaining participants, within stage, into two halves 20
  times (base seed 2041). The halves hold about 710, 650, 500 and 375 participants.
- **Spaces.** `pls_cv` (designated: the case study's configuration, cross-validated on the unit, rank floor 2), the
  predictive subspace of the same model, the fit forced at ranks 2, 8 and 19, and a PCA of the sample at the unit's CV
  rank. Every fit is one `analyze_trajectories` call with `permutations=0`, so no RRPP draw ran.
- **Estimate.** Per point, space and reading, the sum over units of the numerator over the sum of the denominator,
  never clipped (design D4). Readings: the total capture, the net-displacement and interior contrasts, the capture of
  the halves' difference (pure estimation noise), and item 22's two naive captures.
- **Profile `calibration`.** 24 Gaussian cohorts per strength from the full-cohort fit `features ~ stage + covariates`
  and its residual covariance, carrying the cohort's own stage configuration at strengths 0, 0.67 and 1, with 25
  units per cohort and point through the cohort profile's pipeline. Strength 0.67 extrapolates item 22's split-half
  reliability; the pre-flight measured the full cohort's own at 0.663.
- **Code and hosts.** Revision `5cf3b20` from a clean detached worktree on ing, Python 3.11.16, numpy 2.3.5,
  scipy 1.16.3, scikit-learn 1.8.0, pandas 2.3.3, BLAS pinned to one thread. Cohort 32 workers on n3, calibration 64
  on n2, `decide` after both. Launch commands and job times are in `PROVENANCE.txt`.

## 2. Accounting, checks, pre-flights and the reproductions

| profile | units | failed | inconsistent | recorded core-hours | unit median |
|---|---|---|---|---|---|
| cohort | 2,000 / 2,000 | 0 | 0 | 5.77 | 7.1 / 8.4 / 11.6 / 14.5 s at 80 / 200 / 500 / 750 |
| calibration | 7,200 / 7,200 | 0 | 0 | 35.20 | 16.7 s |

- **Every per-unit check passed** on every unit: the anchor's cell table, disjoint and covering parts, orthonormal
  bases at their dimension, Gram traces equal to their scalars, and in the calibration the planted configuration
  recovered from the noise-free features.
- **Pre-flights.** Cohort: replicates 0 and 1 of the full cohort selected rank 8 at mean held-out AUROC 0.599514, and
  the call forced at rank 8 gave the same statistics. Calibration: with the noise scaled to zero, every space's
  estimated capture equals its true capture within 1e-15.
- **Reproductions.** The CV-rank counts at 80, 200 and 500 equal the anchor run's point table, and the median naive
  captures at 750 equal the gain run's cell table (0.53603 and 0.925475 in the fitted span, 0.108803 and 0.548427 in
  the predictive subspace).
- **Cost.** 41.0 recorded core-hours against the smoke's projection of 20.5: units took 1.3 to 2.2 times as long on
  ing as on the workstation.

## 3. The decision

Rule version 1, fixed before any unit ran (design D8); not indicative (500 units per cohort point, 24 cohorts and 25
units in the calibration).

| step | condition | result |
|---|---|---|
| gate: units | failed + inconsistent ≤ 1 % per point | holds (0 everywhere) |
| gate: reproductions | both hold | holds |
| gate: bias | \|bias\| of `pls_cv` total ≤ 0.05 at strengths 0.67 and 1, every point | **fails** at 5 of 8 cells |
| program present | `n80` energy > 2 × its SD at strength 0 | not present: 0.779 against 1.015 |
| resolved | `pls_cv` total half-width ≤ 0.15 | no point: 2.79 / 3.59 / 2.22 / 1.13 |

**Verdict: `gate_failed`.** Under it nothing is read: no capture is stated, no contrast is graded and no follow-up
with plants is licensed. The steps after the gate are recorded in the decision file and quoted here only to show
that the gate was not the only obstacle.

## 4. Why the gate failed

`pls_cv`'s total capture in the calibration, mean over 24 cohorts. The standard error is the cohort SD over √24;
the last column is the chance that an exactly unbiased estimator with that standard error lands within ±0.05.

| strength, point | truth | estimate | bias | SD across cohorts | SE of the mean | bias / SE | P(pass \| unbiased) |
|---|---|---|---|---|---|---|---|
| 1, 80 | 0.321 | 0.244 | **−0.077** | 0.152 | 0.031 | −2.5 | 0.89 |
| 1, 200 | 0.493 | 0.431 | **−0.062** | 0.167 | 0.034 | −1.8 | 0.86 |
| 1, 500 | 0.735 | 0.700 | −0.035 | 0.114 | 0.023 | −1.5 | 0.97 |
| 1, 750 | 0.821 | 0.786 | −0.035 | 0.139 | 0.028 | −1.2 | 0.92 |
| 0.67, 80 | 0.280 | −0.082 | **−0.362** | 1.393 | 0.284 | −1.3 | 0.14 |
| 0.67, 200 | 0.385 | −0.065 | **−0.450** | 1.794 | 0.366 | −1.2 | 0.11 |
| 0.67, 500 | 0.571 | 0.681 | **+0.110** | 1.108 | 0.226 | +0.5 | 0.18 |
| 0.67, 750 | 0.680 | 0.693 | +0.013 | 0.566 | 0.116 | +0.1 | 0.34 |

- **At strength 1 the bias is real and small.** It is negative at every point, about 1.2 to 2.5 standard errors, and
  the PCA and the forced ranks share the sign at every point (−0.067 to −0.110 and −0.030 to −0.086); the predictive
  subspace does not (−0.034 to +0.030).
  Numerator and denominator are unbiased given the fit, so a ratio of sums whose denominator is noisy is the expected
  source (design, Risks); that is a reading for the diagnosis, not a finding of this run.
- **At strength 0.67 the bias cannot be measured to ±0.05.** The cohort-level SD is 0.57 to 1.79, so the 24-cohort
  mean has a standard error of 0.12 to 0.37, and the readings there lie within 1.3 standard errors of zero. Under
  exact unbiasedness the four points would all pass with probability 0.14 × 0.11 × 0.18 × 0.34 ≈ 0.001 (normal
  approximation). The rule fixed the limit before the spread was known, and the spread made the limit unattainable at
  this calibration size.
- **ADNI's program is weaker than strength 0.67.** The cohort's own cross-product energy is 0.70 to 0.78 per point,
  against 0.92 to 0.98 in the calibration at 0.67 and 1.90 to 2.01 at 1: about strength 0.6 in energy. The cohort sits
  where the estimator is noisiest.

## 5. Why it would not have resolved

- **The spread is a cohort's, and one cohort is all ADNI has.** At 80 participants every unit's halves come from a
  remainder of 1,421 of the 1,501 participants. Averaging 500 units and 20 splits removes the noise of the fit and of
  the split, but not that of the remainder's own stage means, which is one draw per cohort. That is the spread the
  calibration measures (design D7), and the interval uses it.
- **The intervals.** `pls_cv`'s total capture reads 0.119, 0.221, 0.366 and 0.477 at the four points, each ± 2.79,
  3.59, 2.22 and 1.13. Every space is unresolved at every point: the narrowest interval is the predictive subspace's
  at 500 participants, ±0.52.
- **Even at strength 1** twice the cohort SD is 0.23 to 0.33, above the limit 0.15.
- **Program presence.** At strength 0 the total energy's SD across cohorts is 0.47 to 0.53, so presence needs about
  1.0. The cohort reads 0.78 at 80 participants; at strength 0.67 the calibration's own cohorts exceed 1.015 about
  half the time (normal approximation). Neither contrast is present at any point: the net displacement's energy is
  0.38 to 0.45 against a null SD of 0.47 to 0.50, the interior's 0.32 to 0.33 against 0.23 to 0.27.

## 6. Descriptive readings (not read under the rule)

These are recorded so that the diagnosis starts from them. The verdict licenses none of them, and none is a capture of
ADNI's program.

**The cohort, by space** (ratio of sums over 500 units; total capture / noise capture / naive capture of the whole
remainder):

| space | 80 | 200 | 500 | 750 |
|---|---|---|---|---|
| `pls_cv` | 0.12 / 0.36 / 0.27 | 0.22 / 0.41 / 0.34 | 0.37 / 0.50 / 0.46 | 0.48 / 0.54 / 0.52 |
| predictive | 0.09 / 0.19 / 0.15 | 0.13 / 0.16 / 0.15 | 0.20 / 0.12 / 0.14 | 0.28 / 0.09 / 0.13 |
| `pca` at the CV rank | 0.13 / 0.51 / 0.37 | 0.18 / 0.55 / 0.42 | 0.22 / 0.60 / 0.50 | 0.28 / 0.62 / 0.55 |
| rank 2 | 0.06 / 0.27 / 0.19 | 0.09 / 0.27 / 0.20 | 0.06 / 0.28 / 0.22 | 0.02 / 0.28 / 0.22 |
| rank 8 | 0.32 / 0.62 / 0.51 | 0.40 / 0.64 / 0.55 | 0.51 / 0.65 / 0.61 | 0.58 / 0.65 / 0.63 |
| rank 19 | 0.61 / 0.76 / 0.70 | 0.69 / 0.78 / 0.75 | 0.75 / 0.79 / 0.78 | 0.79 / 0.79 / 0.79 |

- **The per-unit spread is wide.** For `pls_cv`'s total at 750 participants the unit quantiles run from −3.11 (5 %)
  to 3.42 (95 %), median 0.65.
- **The program is small beside the noise.** The halves' cross-product energy is 0.70 to 0.78, their difference's 5.3
  to 9.9, and the mean split-half correlation 0.25 at 80 participants falling to 0.14 at 750. The pooled cross-Gram
  eigenvalues are 0.63 and 0.15 at 80 participants and 0.58 and 0.12 at 750; the strength-0 calibration gives no
  per-eigenvalue null, so whether the second is above noise is unread.
- **CV ranks** (medians 2, 3, 5, 5) and held-out AUROC (0.52 to 0.58) match the anchor's units, as the reproduction
  requires.

**The calibration's naive capture.** In the Gaussian cohorts the naive capture of `pls_cv` lies within 0.06 of the
truth at every point and both strengths (0.317 / 0.404 / 0.556 / 0.644 against 0.280 / 0.385 / 0.571 / 0.680 at 0.67)
with a cohort SD of 0.08 to 0.11, because there the span holds noise about as much as program (noise capture 0.35 to
0.64 at 0.67). On ADNI that closeness is what was to be measured, so this licenses no reading of item 22's 0.536.

## 7. What follows

- **A diagnosis of the estimator comes before any reading** (the rule's consequence). It needs no new unit: the
  calibration's per-unit records are kept on ing.
    - Split the calibration's spread at strength 0.67 into its between-cohort and within-cohort parts. The interval
      used the SD of a 25-unit estimate; if much of it is within cohorts, the 500-unit cohort estimate is more precise
      than its interval says. If almost all of it is between cohorts, ADNI's one cohort cannot resolve the capture at
      ±0.15 with this estimator, and the number waits for a larger cohort.
    - Locate the bias at strength 1, for example by a delta-method correction of the ratio of sums on the
      calibration's records, and check it against the truth.
    - Any rule version 2 is a new change, written before its units run. A bias limit must be one its calibration can
      measure: at the observed SDs, a standard error of 0.025 at strength 0.67 would need 3,000 to 5,000 cohorts at
      80 and 200 participants.
- **Until then.** The latent plasmode stays open on ADNI. Item 22's held-out readings stay mixtures, the simulated
  0.10 stays the only capture of a weak stage program, and open item 4's candidates cannot use this estimator as a
  pre-check. ROSMAP has no export.

## 8. Limitations

- **The calibration is Gaussian** and its program is the cohort's estimate, about half of it noise in direction; its
  true captures are not ADNI's (design, Risks). The bias it measured belongs to the estimator at ADNI's design, not to
  ADNI's program.
- **Strength 0.67 was a setting.** The cohort's own energy puts it nearer 0.6, below the calibration's interval
  strength.
- **One cohort.** Every descriptive reading in §6 is one cohort's, with dependent units.
- **Section 4's probabilities** use a normal approximation of the 24-cohort mean with the observed SDs; they describe
  the rule's attainability and are not part of the decision.
