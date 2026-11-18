# Notes

## Probe before the proposal (2026-11-18; not a result)

A scratch script outside the repository, on the ADNI export read in place, with the anchor driver's unit builder under a scratch base seed (990118, used by no driver). Nothing was written to the repository or the export. Per point: 32 units, one production fit each (`analyze_trajectories` at the case study's configuration with `permutations=0`), and 20 within-stage splits of the participants outside the sample into two halves. Coordinates: the fitted model's inputs.

The estimate is the pooled ratio over units of `sum(m1·Q·Q'·m2')` to `sum(m1·m2')`, with `m1`, `m2` the halves' centred covariate-adjusted stage configurations and `Q` a space's orthonormal basis. The naive capture is `‖m·Q‖² / ‖m‖²` on the whole remainder's configuration, item 22's reading.

| point | CV rank (median) | fit, s | split-half correlation | cross energy / naive energy | own sample, naive | fitted span: cross / naive | predictive subspace: cross / naive | PCA at the CV rank: cross / naive |
|---|---|---|---|---|---|---|---|---|
| `n80` | 2 | 5.2 | 0.250 | 0.753 / 2.053 | 0.760 | 0.118 / 0.243 | 0.102 / 0.149 | 0.116 / 0.347 |
| `n200` | 3 | 6.2 | 0.238 | 0.782 / 2.211 | 0.874 | 0.224 / 0.355 | 0.063 / 0.104 | 0.237 / 0.398 |
| `n500` | 5 | 8.3 | 0.157 | 0.587 / 2.443 | 0.892 | 0.288 / 0.429 | 0.224 / 0.108 | 0.111 / 0.494 |
| `n750` | 6 | 10.3 | 0.142 | 0.712 / 3.133 | 0.926 | 0.433 / 0.518 | 0.211 / 0.118 | 0.286 / 0.531 |

By stage contrast, in the fitted span (energies with unit-norm contrasts; capture = in-span / total):

| point | net displacement (AD − CN): energy, capture | interior (MCI off the chord): energy, capture | pooled cross-Gram eigenvalues |
|---|---|---|---|
| `n80` | 0.448, 0.041 | 0.305, 0.232 | 0.608, 0.145 |
| `n200` | 0.450, 0.094 | 0.332, 0.399 | 0.616, 0.166 |
| `n500` | 0.278, −0.005 | 0.308, 0.552 | 0.502, 0.085 |
| `n750` | 0.350, 0.354 | 0.362, 0.510 | 0.507, 0.205 |

What the probe shaped:

- **The estimator runs on existing code.** `permutations=0` fits the space and skips RRPP; a unit costs 5 to 10 s at two BLAS threads.
- **A per-unit ratio is unusable.** Its denominator is near zero in some units (per-unit captures from −12.9 to 22.8 at `n500`). Only the ratio of sums over units is read (design D4).
- **The naive capture overstates the cross-product one at every point** (0.518 against 0.433 at `n750`; item 22 read 0.536 there). The remainder's configuration is mostly estimation noise (cross energy 0.59 to 0.78 of a naive 2.05 to 3.13), and that noise lies in the fitted span about as much as in a PCA of the same rank. This motivated the noise-capture reading (design D5).
- **The capture appears to rise with the sample**, from about 0.12 at 80 participants to about 0.43 at 750. The simulated weak design read 0.10.
- **The program may be two-dimensional.** The second pooled eigenvalue is positive at every point, and the two contrasts hold similar energy. Item 22's per-unit eigenvalue ratio (median 0.000) had suggested one dimension. Whether the second is above noise needs the strength-0 calibration (design D7).
- **The net displacement may be captured less than the interior contrast.** At 32 units sharing one cohort this is a lead only.
- **Thirty-two units do not resolve a point.** The points' cross energies differ by 0.2 on one cohort. The run needs the anchor's 500 units per point and a cohort-level standard error, which only a calibration with independent cohorts gives (design D7, D8).

## Smoke (task 6.1, 2026-11-18)

Code: `proposal/adni-latent-plasmode-capture` at `f9d9971` (main) plus this change's uncommitted driver and tests (sha256 prefixes `311f7a8d352ceef9` for `scripts/adni_latent_plasmode_check.py`, `dc9f57a9e1caf0f7` for its test module). Host `nixos-office`, the export read in place (`../adni_sexdiff/results/motco_input`). Python 3.11.15, numpy 2.3.5, scipy 1.16.3, scikit-learn 1.8.0, pandas 2.3.3, motco 0.25.0 from the checkout's `src`. BLAS pinned to one thread in the shell; eight workers. Output off the repository: `~/motco-runs/adni-latent-plasmode-smoke`. The docstring's seven commands, run as written (`--reps 4`; `--cohorts 2 --units 2`).

- **Cohort pre-flight passed** (34 s). On cohort replicates 0 and 1 the cross-validated fit selected rank 8 at a mean held-out AUROC of 0.599514, and the call forced at rank 8 gave the same statistics (CV 16.1 s and 16.2 s; forced 0.03 s). No permutation was drawn and no statistic is recorded.
- **The full cohort's reliability** (recorded, not judged; 20 within-stage splits into halves of 750 and 751): mean split-half correlation 0.282 (range −0.057 to 0.568), Spearman–Brown 0.440 at the full cohort, square root 0.663 — beside the calibration's preset strength 0.67 (from item 22's 0.166 at 375 per half). Pooled cross-Gram eigenvalues 0.591, 0.220, 0.
- **Calibration pre-flight passed** (27 s). At strength 1 with the noise scaled to zero, one unit per point: every space's estimated capture equals its true capture (total and both contrasts) within 9e-16; the CV selected rank 6 at every point. sklearn's "y residual is constant" warnings on the noise-free cohort are silenced in the driver's fits.
- **Runs.** Cohort: 16 of 16 units (four per point) completed, 0 failed, 0 inconsistent, no subsample or label redraw; 22 s wall. Calibration: 48 of 48 (two cohorts × two units × three strengths × four points), 0 failed, 0 inconsistent; 55 s wall. Every check passed on every unit, `truth` included (deviations about 1e-13).
- **The units are the parents' units.** The four `n750` units reproduce the gain run's per-unit records (`ing:/home1/tgonza/MOTCO-adni-gain-run/out/plants/merged.jsonl`, anchor cell, replicates 0–3, read only): CV ranks 3, 5, 7, 6, and all four naive captures (remainder and sample, fitted span and predictive subspace) to within 1e-15 (the last digits differ between the two machines). The merge's reproduction against the committed 500-unit tables fails at four units, as it must; the smoke's decision is `gate_failed` for that and for the two-cohort calibration's bias.
- **Records.** No record holds a list of a sample's, remainder's, half's, the cohort's or the feature set's length; the export's five files hash as the settings bind them after the run. (ADNI's identifiers are integers, so a text search of the aggregates is uninformative; the structural check is the one that holds.)
- **Timing** (median seconds per unit, eight workers): cohort `n80` 5.6, `n200` 6.5, `n500` 8.7, `n750` 10.8; calibration 5.6, 6.6, 9.0, 11.0 (the same at each strength). Almost all of it is the one cross-validated fit; the three forced fits take about 0.05 s together and the configurations and splits 0.1 s.
- **Projection** (the committed counts at those medians): cohort 500 × (5.6 + 6.5 + 8.7 + 10.8) s = **4.4 core-hours**; calibration 1,800 × (5.6 + 6.6 + 9.0 + 11.0) s = **16.1 core-hours**; **20.5 recorded core-hours in all** (the proposal estimated about 25), about 2.6 h at eight workers or 1.9 h at eleven on this machine (12 physical cores), plus two minutes of pre-flights.
- **Descriptive readings at four units** (none decides anything; the smoke's intervals rest on two cohorts): the pooled `pls_cv` total capture read 0.22 at `n80`, 0.05 at `n200`, 0.57 at `n500` and 3.4 at `n750`, where the four units' mean cross energy was −0.07. At `n750` the halves (375 each) leave the denominator dominated by estimation noise unit by unit, as the probe found; whether 500 units and the calibration's spread resolve it is the run's question. The halves' difference (`noise`) sat 0.44 to 0.56 inside the fitted span at every point.
