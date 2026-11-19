# Proposal

## Why

How much of a real stage program the production span sees is not established (roadmap open item 2, the latent plasmode; `docs/roadmap.md`, "Limits of the span measurement"). The only capture of a weak stage program is simulated: 0.10 at SEA-AD's stage signal (item 14).

Item 22 came closest. It read the share of a held-out stage configuration inside the span fitted on 750 ADNI participants: 0.536. That number is not the answer, because the held-out configuration is about 0.3 stage program and 0.7 estimation noise (Spearman–Brown 0.285), and a share of the mixture cannot be split into its two parts. The same noise made the held-out plants unreadable for specificity.

Two independent halves of the held-out participants remove the noise in expectation. With `m1` and `m2` their centred stage configurations and `Q` the span fitted on a third, disjoint sample, `tr(m1·Q·Q'·m2')` and `tr(m1·m2')` are unbiased for the true program's energy inside the span and in total. ADNI's cells can be split three ways (the smallest holds 125 participants); SEA-AD's cannot (9).

A probe on scratch seeds (2026-11-18, `notes.md`; not a result) found the estimator feasible and the naive reading biased upward: about 0.43 against a naive 0.52 at 750 participants, and about 0.12 at 80. It also found that 32 units do not resolve a point and that a per-unit ratio is unusable.

## What Changes

- **A driver, and no package change.** `scripts/adni_latent_plasmode_check.py` (`preflight` | `run` | `merge` | `decide`, with `--profile cohort | calibration`). Every fit is one `analyze_trajectories` call on existing options with no permutations. `motco` stays at 0.25.0.
- **Profile `cohort`: the capture on ADNI's blocks.** It imports the anchor driver's contract, subsampling and seeds unchanged, at 80, 200, 500 and 750 participants, 500 units each. The units at 80, 200 and 500 are the anchor run's and those at 750 the gain run's.
  - Each unit fits the spaces on its sample and splits the participants outside it, within stage, into two halves, 20 times.
  - The capture of the covariate-adjusted pooled stage program is the ratio, summed over units, of the halves' cross-product inside a space to their cross-product in total.
  - It is read in total and by stage contrast: the net displacement, and the interior stages' deviation from it.
  - The capture of pure estimation noise (the halves' difference) is read beside it.
- **Spaces.** Production's fitted span at the CV rank; the predictive subspace of the same fit; the fit forced at ranks 2, 8 and 19; and a PCA of the sample at the CV rank as the label-free reference.
- **Profile `calibration`: the estimator with the truth known.** Gaussian cohorts at ADNI's design, covariates and residual covariance, carrying the cohort's own stage configuration at strengths 0, 0.67 and 1. It gates the estimator's bias, gives the cohort-level standard error of the estimate, and gives the level a contrast's energy must exceed to count as present.
- **A decision fixed before any unit runs.** `gate_failed`, `no_program`, `unresolved`, `partly_measured` or `measured`, with a grade for each stage contrast and a statement of what each licenses.
- **No statistic of group difference is computed.** No RRPP draw, no plant, and no analysis under the real sex labels. The stage program is pooled over groups, and the fit is group-blind.
- **Report and records.** A dated report, the aggregates under `results/`, the reading on the measurement-visibility page and the verdict in the roadmap.

Not in this change: plants along a denoised program and their power (a follow-up the decision licenses or not); the split test's reading on ADNI's out-of-span cells; open item 4's candidate spaces, which can import this driver's estimator; ROSMAP.

## Capabilities

### New Capabilities

- `cohort-stage-program-capture`: the share of a cohort's real stage program inside a fitted measurement space, estimated from independent held-out halves on the cohort's exported blocks, calibrated on synthetic cohorts with the truth known, and read under a predeclared decision.

### Modified Capabilities

None. `latent-space-metric-compatibility` owns the measurement-visibility page; its requirements do not change, and the page gains a reading as it did for items 21 and 22 through a task.

## Impact

- **Code:** `scripts/adni_latent_plasmode_check.py` and `tests/test_adni_latent_plasmode_check.py`. It imports `adni_anchor_check` and `adni_gain_specificity_check` unchanged and binds their committed tables by sha256.
- **Package, API, CLI, dependencies:** none. No version bump and no changelog entry.
- **Data:** the ADNI export is read in place on a machine that holds it and never copied. Records hold seeds and statistics only.
- **Compute:** about 25 core-hours on one machine by the probe's timing (2,000 cohort units and 7,200 calibration units at 5 to 10 s each). The smoke measures it, and the full runs wait on the user's approval of the cap.
- **Docs:** `docs/reports/adni-latent-plasmode-<date>.md`, `docs/measurement-visibility.md`, `docs/roadmap.md`, `mkdocs.yml`, and one row in `AGENTS.md` (Validation drivers).
