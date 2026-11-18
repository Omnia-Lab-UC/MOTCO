# Notes

## Exploratory pass before the proposal (2026-11-18; not a result)

A scratch script outside the repository, run on ing (srun jobs 945689 and 945690) on the latent-plasmode run's merged records. It used the parent driver's `unit_reading` and `truth_reading` and wrote nothing to the repository or the run's output. The rule in `design.md` was written after it.

What it found, on `pls_cv`'s total capture:

- **The parts are unbiased.** Means of `N − TN` and `D − TD` lie within 0.7 cohort-level standard errors of zero at every gate cell.
- **The aggregate is on its truth.** Over all 600 units of a cell, the ratio of sums minus the truth's has a 95 % cluster-bootstrap interval covering zero at every gate cell (strength 1: −0.028 to +0.003).
- **The spread is between cohorts.** The between share of the variance of a cohort's mean `N` and `D` is 0.89 to 0.99. The delta-method bias at 500 units is about the same as at 25 (−0.024 to −0.049 at strength 1).
- **At strength 0.67 the per-cohort SD is driven by one or two cohorts** with near-zero denominators. Robust SD 0.17 to 0.30, median error −0.03 to −0.05. A fit of `var(m) = b + w/m` on subsets of units was unstable there (one subset SD of 14,858) and is not carried into the driver.
- **ADNI's within-unit SDs of `N` and `D`** are 10 to 28 % below the calibration's at strength 0.67, and close to strength 0's.
- **ADNI's `t = D̄ / √V_DD`** is 1.22 to 1.55 under the raw stand-ins and 1.52 to 1.90 under the scaled ones, all below 1.96. So no 95 % Fieller interval is bounded at any point.
