# The ADNI latent plasmode, diagnosed — report (2026-11-18)

**Change:** `adni-latent-plasmode-diagnosis` (MOTCO 0.25.0, no package change). **Reads:**
`results/adni-latent-plasmode-diagnosis-2026-11-18/` (`adni_latent_plasmode_diagnosis.json`, `component_table.csv`,
`spread_table.csv`, `adni_table.csv`, `reproduction.json`, `settings.json`, `PROVENANCE.txt`). **Parent:** the
[ADNI latent plasmode](adni-latent-plasmode-2026-11-18.md) (its committed settings, tables and decision, bound by
sha256). **Concept page:**
[measurement visibility](../measurement-visibility.md#the-latent-plasmode-on-adni-gate-failed).

The latent plasmode estimated the share of ADNI's covariate-adjusted stage program inside production's fitted span
from two independent held-out halves. It read `gate_failed`: in its calibration on synthetic cohorts, the bias of
`pls_cv`'s total capture exceeded ±0.05 at five of eight strength × point cells. Its rule names the consequence: a
diagnosis of the estimator comes before any reading. This diagnosis re-reads the run's 9,200 records (no new unit) and
asks where the bias comes from, where the spread lives, and whether any interval for ADNI's capture is bounded.

**In one paragraph.** **The reading is `information_limited`.**

- **The estimator is sound.** A unit's numerator and denominator are unbiased for their truths, and pooled over a
  cell's 24 cohorts the ratio sits on its truth. The gate's bias is the nonlinearity of a ratio with a noisy
  denominator, within one cohort.
- **The spread belongs to the cohort.** 89 to 99 % of the variance of a cohort's mean denominator is between cohorts,
  so more units or splits of one cohort do not narrow it.
- **ADNI's program does not clear that spread.** Its mean denominator is 1.22 to 1.90 of its standard deviations from
  zero, under every stand-in for the cohort-level variance. So no two-sided 95 % interval for the capture is bounded at
  any point.
- **What is left: one-sided bounds** (recorded, not judged). A capture lies in [0, 1], and within that range ADNI's
  sets give upper bounds: production's span holds at most 0.46 to 0.57 of the program at 80 participants and 0.52 to
  0.65 at 200, depending on the stand-in.
- **The consequence.** No capture is stated on ADNI, and the parent's `gate_failed` stands. A resolved capture waits
  for a larger cohort or several cohorts.

**How to read the verdict.** The rule (version 1) was written after an exploratory pass on the same records (the
change's `openspec/changes/archive/2026-11-18-adni-latent-plasmode-diagnosis/notes.md`). It is a reading with its thresholds stated, not a predeclared test. The range restriction was added after the first
recorded run, and it is not judged.

## 1. Configuration and provenance

- **Inputs.** The run's merged records on ing, read in place: 7,200 calibration and 2,000 cohort records (sha256 in
  `reproduction.json`). Every record carries its profile's committed settings hash, and each expected unit is present
  once. No export is read.
- **Reproduction.** The parent driver's own table builders, run on the records, write the committed
  `calibration_table.csv` and `cohort/capture_table.csv` byte for byte.
- **Code.** Revision `e9ec526` from a clean worktree on ing, srun job 945803, about four minutes on two cores. Two
  earlier recorded runs are in `PROVENANCE.txt` and the change's notes. They had the same tables; one dropped the
  strength-0 stand-ins (a bug, fixed with a test), and neither carried the ranges.
- **The designated reading** is `pls_cv`'s total capture. The tables carry every space (the predictive subspace, a PCA
  at the CV rank, forced ranks 2, 8 and 19) and both stage contrasts.

## 2. The reading

| finding | holds when (`pls_cv`, total) | result |
|---|---|---|
| F1, the parts | \|z\| of `N − TN` and `D − TD` below 3 at every gate cell | holds: largest \|z\| 0.67 |
| F2, the aggregate | \|z\| of the pooled ratio's bias below 3 at every gate cell | holds: largest \|z\| 1.10 |
| F3, the spread | between share of `var(D̄)` at strength 0.67 at least 0.5 at every point | holds: 0.89 to 0.99 |
| F4, ADNI's interval | unbounded at every point under the primary stand-in (0.67, raw) | holds; all four stand-ins agree |

The verdict order is `reproduction_failed`, `estimator_biased` (F1 or F2), `resolvable` (F4), `units_limited` (F3),
then `information_limited`. Across every space and reading at the gate strengths, the largest \|z\| of either part or
of the aggregate is 1.75.

## 3. Where the bias comes from (F1, F2)

`pls_cv`'s total capture in the calibration. The parts' bias has a standard error from the 24 cohort means. The pooled
ratio is over all 600 units of a cell, with a cluster-bootstrap standard error and 95 % interval over cohorts.

| strength, point | bias of `N` (z) | bias of `D` (z) | pooled ratio − truth | its 95 % interval | the gate's per-cohort bias |
|---|---|---|---|---|---|
| 0.67, 80 | +0.049 (0.53) | +0.063 (0.54) | +0.033 | −0.101 to +0.143 | −0.362 |
| 0.67, 200 | +0.048 (0.47) | +0.046 (0.40) | +0.032 | −0.097 to +0.144 | −0.450 |
| 0.67, 500 | +0.029 (0.25) | +0.023 (0.20) | +0.017 | −0.109 to +0.111 | +0.110 |
| 0.67, 750 | +0.015 (0.13) | +0.022 (0.18) | +0.000 | −0.095 to +0.073 | +0.013 |
| 1, 80 | −0.064 (−0.67) | −0.021 (−0.13) | −0.028 | −0.084 to +0.016 | −0.077 |
| 1, 200 | −0.024 (−0.19) | −0.018 (−0.11) | −0.008 | −0.066 to +0.042 | −0.062 |
| 1, 500 | −0.015 (−0.10) | −0.029 (−0.18) | +0.003 | −0.038 to +0.037 | −0.035 |
| 1, 750 | −0.065 (−0.44) | −0.087 (−0.60) | +0.003 | −0.029 to +0.029 | −0.035 |

- **The parts are unbiased**, as the design required: the halves are independent of each other and of the fit.
- **Pooled over cohorts, the ratio is on its truth.** The gate judged the mean over cohorts of each cohort's own ratio.
  With 25 units of one cohort that ratio's denominator stays noisy, and a ratio's error does not average out when its
  denominator is noisy.
- **At strength 1 the delta method accounts for the size.** Its first-order bias at 25 units per cohort is −0.026 to
  −0.050, against the observed −0.035 to −0.077. At strength 0.67 the delta method does not apply: the cohort mean
  denominator's coefficient of variation is 0.59 to 0.66.
- **At strength 0.67 one or two cohorts dominate the gate's numbers.** One or two of the 24 cohorts have a mean
  denominator below a quarter of the cell's. They inflate the per-cohort SD to 0.57 to 1.79, while the robust SD
  (IQR / 1.349) is 0.20 to 0.31 and the median error −0.026 to −0.053.

## 4. Where the spread lives (F3)

| strength, point | between share of `var(N̄)` | between share of `var(D̄)` | delta-method bias at 25 units | at 500 units |
|---|---|---|---|---|
| 0.67, 80 | 0.981 | 0.991 | −0.147 | −0.147 |
| 0.67, 200 | 0.979 | 0.981 | −0.146 | −0.145 |
| 0.67, 500 | 0.959 | 0.952 | −0.160 | −0.156 |
| 0.67, 750 | 0.899 | 0.893 | −0.118 | −0.109 |
| 1, 80 | 0.968 | 0.994 | −0.044 | −0.044 |
| 1, 200 | 0.971 | 0.989 | −0.050 | −0.049 |
| 1, 500 | 0.963 | 0.966 | −0.034 | −0.033 |
| 1, 750 | 0.915 | 0.912 | −0.026 | −0.024 |

- **More units would not help.** Every unit's halves come from the participants outside its sample, and at 80
  participants that is 1,421 of the same 1,501. The units and splits average out the fit and the split, but not the
  noise of the cohort's own stage means. So ADNI's 500 units are about as precise as the calibration's 25 per cohort,
  and the ratio's bias barely moves from 25 to 500 units.
- **The share is lowest at 750 participants,** where the remainders overlap least.
- **Across every space and reading at strength 0.67** the between share of `var(D̄)` is 0.87 to 0.99.

## 5. ADNI's own interval (F4)

- **ADNI's inputs.** For each point, ADNI's mean numerator and denominator over its 500 units.
- **Their covariance.** `V = f·B + W_ADNI / 500`. ADNI's between-cohort covariance cannot be observed from one cohort,
  so `B` comes from the calibration at strength 0 or 0.67. It is used raw (`f = 1`) or scaled by ADNI's within-cohort
  variance of the denominator over the calibration's (`f` 0.52 to 0.82).
- **When the interval is bounded.** Fieller's 95 % interval for the ratio is bounded only when `t = D̄ / √V_DD`
  exceeds 1.96.

| point | ADNI's mean `D` | capture (ratio of sums) | `t`: 0.67 raw / 0.67 scaled / 0 raw / 0 scaled | upper bound within [0, 1], same order |
|---|---|---|---|---|
| 80 | 0.779 | 0.119 | 1.36 / 1.90 / 1.54 / 1.75 | 0.57 / 0.49 / 0.49 / 0.46 |
| 200 | 0.769 | 0.221 | 1.38 / 1.85 / 1.51 / 1.67 | 0.65 / 0.59 / 0.54 / 0.52 |
| 500 | 0.710 | 0.366 | 1.26 / 1.52 / 1.55 / 1.70 | 0.77 / 0.74 / 0.73 / 0.69 |
| 750 | 0.699 | 0.477 | 1.22 / 1.58 / 1.41 / 1.70 | 0.81 / 0.77 / none / none |

- **No interval is bounded.** `t` stays below 1.96 at every point under every stand-in. ADNI's program energy (0.70 to
  0.78 in pooled SDs) is about 1.2 to 1.9 of its own cohort-level standard deviations from zero.
- **What the sets do say.** Every set is the line outside two roots, or the whole line. Within the capture's range
  [0, 1], all but two are an interval from 0 to the bound in the last column. At 750 participants under the strength-0
  stand-ins the set covers the whole range.
- **The bounds** are recorded, not judged, and they are bounds, not estimates. At 80 and 200 participants production's
  span holds at most about half to two thirds of ADNI's program. For comparison at 80 participants, the label-free PCA
  is bounded at 0.59 to 0.77 and the predictive subspace at 0.30 to 0.37 (`adni_table.csv`).

## 6. What follows

- **No capture is stated on ADNI.** `gate_failed` stands, and so does the simulated 0.10 as the only capture of a weak
  stage program. A rule version 2 on ADNI that resolves the capture at ±0.15 would fail for the same reason: the
  cohort's own spread, not the units or the gate.
- **The estimator stands.** Its parts and aggregate are unbiased, so it is fit for a cohort whose program energy
  clears its cohort-level spread. ROSMAP is the next candidate once an export exists; several cohorts pooled would also
  shrink the spread.
- **One-sided bounds are the most ADNI can give.** A rule that states them (for example, "production's span holds at
  most about half of the program at 80 participants") would be its own change, written before it reads, with the
  stand-in fixed in advance.
- **Item 22's held-out readings stay mixtures,** and open item 4's candidates have no real-noise capture pre-check on
  ADNI.

## 7. Limitations

- **The reading was written after the numbers were seen,** and the range restriction after the first recorded run. Both
  are recorded as such. The parent's verdict is the only predeclared one, and it is unchanged.
- **The cohort-level covariance is a stand-in** from Gaussian synthetic cohorts. Its 24 cohorts estimate a variance
  with about 30 % relative error, and ADNI's within-cohort variance is 31 to 48 % below the calibration's at strength
  0.67 (`within_ratio_D`). The four stand-ins bracket that range and agree on F4. The bounds in §5 move with them.
- **Fieller's interval assumes the cohort means are close to normal.** The strength-0.67 calibration shows heavy tails
  near a zero denominator.
- **One cohort.** ADNI's between-cohort spread cannot be observed, only modelled.
