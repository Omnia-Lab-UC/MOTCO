# The `angle` RRPP null is strongly non-pivotal — and that is what makes the test work

**Run date:** 2026-09-01
**Code revision:** `899179586c6c04c97e3485bcf19e504309320ba8` (branch `feat/diagnose-angle-null-pivotality`,
parent `bc0d42c`, the merge of the orientation sign-anchor fix)
**Configuration:** `examples/trajectory_power_study/angle_pivotality_diagnostic.json`
**Records:** `results/angle-pivotality-2026-09-01/merged.jsonl` (500 replicates, 0 failures, 199 permutations each)
**Tables:** `results/angle-pivotality-2026-09-01/report/` — `pivotality_association.csv`,
`pivotality_rejection_split.csv`, `pivotality_standardized.csv`
**Provenance:** `results/angle-pivotality-2026-09-01/PROVENANCE.txt`
**Change:** `openspec/changes/diagnose-angle-null-pivotality`
**Answers:** [Phase 5 readiness](https://github.com/Omnia-Lab-UC/MOTCO/blob/c440070dcb3c5bcc6ddd2227cf7adbef18573b6a/docs/phase5-readiness.md) item 1 (blocking)
**Environment:** Python 3.11.15, numpy 2.3.5, pandas 2.3.3, scikit-learn 1.8.0, scipy 1.16.3; no R at runtime

## Summary

The Phase 4 hypothesis is **confirmed**: within a replicate, the RRPP permutation null for `angle` moves
almost one-for-one with that replicate's own observed angle. Regressing each replicate's own 95th-percentile
null on its observed statistic gives a slope of **0.811** in the orientation cell at effect 1.00 and
**0.87–0.96** in every other cell measured. The statistic is not pivotal, and the effect is large.

This fully accounts for the rejection inversion. Among the 100 orientation replicates at effect 1.00, the
32 that **fail** to reject have a *larger* mean observed angle (60.8°) than the 68 that reject (46.5°) — but
their own critical values differ by far more: **103.6° versus 15.8°**. The outcome is decided by the
critical value, not by the observed statistic, and the pattern is confined to `angle`.

**The proposed remedy does not follow.** The readiness worklist anticipated that a confirmed association
would mean "a pivotal statistic or a studentized test." Both are ruled out by measurement:

- **Within-replicate studentization is a no-op.** Standardizing observed and null by the same constants is a
  strictly monotone transform of both sides of the comparison; the permutation p-value is identical. This is
  a theorem, not an empirical finding, and it is pinned by
  `tests/test_pivotality.py::test_within_replicate_studentization_leaves_the_p_value_unchanged`.
- **Cross-replicate standardization gains nothing.** Recalibrating each replicate as
  `z = (observed − null_mean) / null_sd` against the pooled null-control `z` distribution moves orientation
  `angle` power from 0.68 to **0.70** — two replicates, against a Monte Carlo standard error of 0.047.
- **Freezing the null is catastrophic.** Replacing the per-replicate critical value with a single fixed
  threshold calibrated on the null cell (its 95th-percentile observed angle, 153.2°) drops orientation
  `angle` power from 0.68 to **0.01**.

The per-replicate null is not a defect that costs power. It is adaptation to a latent-space noise level that
varies enormously across replicates, and it is worth roughly 67 points of power relative to any fixed
threshold. The tracking is the price of that adaptation.

**Verdict: `angle` is strongly non-pivotal, the non-pivotality explains the inversion, and it is
load-bearing rather than corrigible.** The remaining shortfall is a property of the design point, not of the
statistic or its calibration.

## Provenance: this reads the same records that produced the shortfall

The diagnostic was configured to hold every generator, integration, and evaluation parameter at the Phase 4
pilot's values, and it does so exactly. Comparing per replicate against the sign-fix re-run
(`results/orientation-signfix-2026-08-28/merged.jsonl`):

- **400 of 500 records match on all four cells they share** — identical generator seeds, identical `delta`,
  `angle`, and `shape` statistics to 1e-8 relative, identical p-values to 1e-12.
- Orientation `angle` power reproduces at **0.68**, `delta` at 0.65, `shape` at 0.99 — the sign-fix
  operating point, replicate for replicate.
- The fifth cell (`magnitude` at effect 1.00, absent from the sign-fix re-run) matches **Phase 4** on all 100
  generator seeds and on `delta` and `shape` in 100 of 100 replicates; `angle` differs in 15 of 100, which is
  the sign-anchor correction and nothing else.

So the numbers below are not a re-derivation on fresh draws. They are a decomposition of the same replicates
that yielded 0.68.

### Cells

`enumerate_study` emits five cells, not the three requested: it unconditionally adds a `none` Type I
baseline and a `translation` negative control. Both are retained deliberately — they supply the null-control
reference distribution the cross-replicate counterfactual calibrates against.

| Cell | Phase | n | Role |
|---|---|---|---|
| `orientation` @ 1.00 | `power_primary` | 100 | the shortfall |
| `magnitude` @ 1.00 | `power_primary` | 100 | comparator; its `delta` reaches power 1.00 |
| `translation` @ 1.00 | `power_primary` | 100 | null control that behaved normally in Phase 4 |
| `translation` @ 1.00 | `type_i_baseline` | 100 | mandated negative control |
| `none` | `type_i_baseline` | 100 | mandated zero-effect baseline |

## 1. The association

Slope of each replicate's own null 95th percentile on its own observed statistic. Both axes carry the same
units within a statistic, so the slope reads directly: *how much of a one-unit increase in the observed
statistic is matched by an increase in the critical value it must beat.*

| Cell | `delta` | `angle` | `shape` |
|---|---|---|---|
| `none` (zero effect) | 0.339 | **0.931** | 0.449 |
| `translation` @1.00 (`type_i_baseline`) | 0.237 | **0.925** | 0.373 |
| `translation` @1.00 (`power_primary`) | 0.152 | **0.960** | 0.348 |
| `magnitude` @1.00 | **0.030** | 0.874 | −0.000 |
| `orientation` @1.00 | 0.018 | **0.811** | **0.058** |

Bold marks each statistic in the cell where it carries signal, plus `angle` throughout.

**Read the columns, not the rows.** Under the zero-effect null every statistic tracks its own null to some
degree, and that is expected: observed and permuted values are exchangeable draws from the same
distribution, so a replicate with an unstable latent space produces both a large observed statistic and a
large null. Tracking under the null is not evidence of a problem.

What separates a usable statistic is what happens when signal is added. For `delta`, the slope in its own
cell **collapses** from 0.339 to 0.030; for `shape`, from 0.449 to 0.058. The signal moves the observed
statistic while leaving the null where it was — which is exactly what a test needs. For `angle` in its own
cell the slope barely moves, 0.931 → 0.811. Roughly four fifths of any orientation signal that reaches the
observed angle also reaches the bar it has to clear.

Correlations, with Fisher-z 95% intervals (`pivotality_association.csv`, `null_target = q95`), confirm the
association is not a small-sample artifact — every interval excludes zero:

| Cell | `angle` correlation | 95% interval |
|---|---|---|
| `none` | +0.634 | [+0.500, +0.738] |
| `translation` @1.00 (`type_i_baseline`) | +0.609 | [+0.469, +0.719] |
| `translation` @1.00 (`power_primary`) | +0.562 | [+0.411, +0.683] |
| `orientation` @1.00 | +0.526 | [+0.367, +0.655] |
| `magnitude` @1.00 | +0.473 | [+0.305, +0.613] |

The same table reports the association against the null's **mean** and **sd**, which track more weakly
(orientation `angle`: slopes 0.306 and 0.349, correlations +0.585 and +0.563). The upper tail moves more
than the centre — which is why the summary retains quantiles and not just moments: the alpha-level critical
value is not recoverable from the mean and sd of a bounded, skewed angular null.

## 2. It accounts for the inversion

From `pivotality_rejection_split.csv`, orientation at effect 1.00, α = 0.05:

| Statistic | n reject | n fail | mean observed (reject) | mean observed (fail) | mean critical (reject) | mean critical (fail) |
|---|---|---|---|---|---|---|
| `angle` | 68 | 32 | 46.5° | **60.8°** | 15.8° | **103.6°** |
| `delta` | 65 | 35 | 3.069 | 0.547 | 1.025 | 1.048 |
| `shape` | 99 | 1 | 0.0705 | 0.0103 | 0.0233 | 0.0266 |

The `angle` row is inverted — non-rejecting replicates carry the larger observed statistic — and it is the
only inverted row in any signal-carrying cell. For `delta` and `shape` the critical value is essentially
constant across the split (ratios 1.02× and 1.14×) while the observed statistic separates by 5.6× and 6.9×.
For `angle` the observed statistic separates by only 1.3× while the critical value separates by **6.6×**.

Per-replicate `angle` critical values in this cell span **5.0° to 176.6°**, a 35-fold range. That dispersion,
not the observed angle, decides the outcome.

The Phase 4 report recorded this inversion pre-sign-fix as 66.7° versus 52.5°. Post-fix it is 60.8° versus
46.5°: the inversion survived the correction, exactly as the sign-anchor report predicted, and it now has a
mechanism.

One further inverted row appears — `angle` in the `translation` Type I baseline (6.9° rejecting vs 14.0°
non-rejecting). It rests on 4 rejections out of 100 and should be read as noise, not as a second finding.

## 3. Why no remedy recovers the power

### Within-replicate studentization is a no-op

Dividing the observed statistic and every null draw by the same constant is strictly monotone on both sides
of the comparison that defines the permutation p-value. The p-value cannot change. This closes off the
reading that a "studentized angle test" is available but unimplemented.

### Cross-replicate standardization: +0.02

`pivotality_standardized.csv` reports, per cell and statistic, the as-specified rejection rate beside the
rate obtained by comparing each replicate's `z = (observed − null_mean) / null_sd` against the 95th
percentile of the pooled null-control `z` distribution:

| Cell | Statistic | As specified | Standardized | Control? |
|---|---|---|---|---|
| `orientation` @1.00 | `angle` | 0.68 | **0.70** | no |
| `orientation` @1.00 | `delta` | 0.65 | 0.64 | no |
| `orientation` @1.00 | `shape` | 0.99 | 0.99 | no |
| `magnitude` @1.00 | `delta` | 1.00 | 1.00 | no |
| `magnitude` @1.00 | `angle` | 0.00 | 0.00 | no |
| `none` | `angle` | 0.05 | 0.05 | yes |
| `translation` @1.00 (`type_i_baseline`) | `angle` | 0.04 | 0.04 | yes |
| `translation` @1.00 (`power_primary`) | `angle` | 0.03 | 0.06 | yes |

The controls hold their nominal level under both rules, so the counterfactual is calibrated. It simply does
not help: two replicates on orientation `angle`, against a Monte Carlo standard error of 0.047. The reason
is visible in §1 — the null's *sd* tracks the observed statistic too (slope 0.349), so standardizing removes
the nuisance and a comparable share of the signal with it.

> **The cross-replicate counterfactual is a diagnostic, not a deployable test.** It borrows a reference `z`
> distribution from null-control cells, which do not exist alongside a real dataset. It is reported to
> measure how much power the tracking costs, never as a procedure to adopt.

### A fixed threshold destroys the test

The null cell's observed `angle` distribution has a median of **5.1°** and a 95th percentile of **153.2°** —
it is extraordinarily heavy-tailed. Using that 95th percentile as one fixed critical value for every
replicate, orientation at effect 1.00 rejects **1 time in 100**, against 0.68 as specified. Only 1 of the 32
non-rejecting replicates carries an observed angle above it.

This is the finding that reframes the item. The per-replicate null is not a tax on power; it is worth about
67 points of it. A replicate-specific critical value is the only thing that makes a latent angular statistic
testable at all when latent-space stability varies this much across draws.

## 4. What the shortfall actually is

The 32 non-rejecting replicates are the high-noise draws: their latent trajectory geometry is unstable
enough that a real effect-1.00 orientation change (cell median observed angle 41.8°, IQR [23.5°, 70.3°])
falls inside their own permutation null. The test is behaving correctly — it is declining to call an angle
significant when that angle is not resolvable in the space it was measured in.

**The instability is invisible to every integration diagnostic currently persisted.** Within the orientation
cell:

- Selected PLS dimensionality is 3 in 93 of 100 replicates (2 in the remaining 7).
- CV mean AUROC is 1.0000 in **every** replicate, rejecting and non-rejecting alike.
- Correlation of log(null q95) with selected dimensionality is **+0.139**; with CV mean AUROC, **−0.179**.

So the driver of the null width is not latent dimensionality and not stage-separation quality. Nothing the
harness records today identifies which replicates will be resolvable, which means Phase 5 cannot condition
on it and cannot report it as a covariate.

## 5. Phase 5 consequence

> **`angle` proceeds as specified.**

The statistic and its permutation test are sound. Type I level is nominal in every control cell under both
decision rules (0.03–0.05 for `angle`), RRPP p-values remain exact under exchangeability whether or not the
statistic is pivotal, no substitute statistic is available, and the two candidate remedies were measured and
found to be respectively inert and destructive. There is no calibration defect to fix and nothing to
replace.

This decision deliberately does **not** settle the 0.80 power floor. This diagnostic measured one design
point (n = 300, four stages, ~660 features, PLS at 3 latent dimensions) and cannot say whether the gap
closes with more samples per group-stage cell. That question belongs to readiness item 4, and it now has a
concrete lever: the dispersion of the per-replicate `angle` null. Power for orientation is governed by how
tightly latent trajectory geometry is determined by the data, so the design-point study should measure how
the spread of `null_summary["angle"]["q95"]` contracts with samples per group-stage cell, not just how the
rejection rate moves.

### Consequences for the other readiness items

- **Item 2 (orientation → shape).** Unchanged in substance, but the diagnostic removes one candidate
  explanation: `shape` in the orientation cell is nearly pivotal (slope 0.058, critical values 0.0233 vs
  0.0266 across the rejection split), so its 0.99 rejection rate is not a null-tracking artifact. The
  response is real in the latent space and item 2's framing as a projection-versus-construction question
  stands.
- **Item 3 (latent dimensionality).** Cannot be settled with the diagnostics currently recorded. Selected
  dimensionality is effectively constant and CV AUROC is saturated at 1.0, and neither correlates with the
  null width. Re-sizing the latent space needs a direct measure of latent trajectory-geometry stability,
  which the harness does not yet compute.
- **Item 4 (design point).** Inherits the 0.80-floor question, with the null-dispersion lever above.

## Reproducing

```bash
# ~26 minutes on 16 cores; 500 work units. Do NOT pass --n-jobs.
for i in $(seq 0 15); do
  uv run python scripts/run_study_shard.py \
    --config examples/trajectory_power_study/angle_pivotality_diagnostic.json \
    --out-dir results/angle-pivotality-2026-09-01 \
    --shard-index "$i" --n-shards 16 --error-policy record &
done
wait
uv run python scripts/motco_study.py merge --out-dir results/angle-pivotality-2026-09-01
uv run python scripts/angle_null_pivotality.py \
  --merged  results/angle-pivotality-2026-09-01/merged.jsonl \
  --out-dir results/angle-pivotality-2026-09-01/report
```

The shard and merged JSONL are gitignored as regenerable; the `report/` tables and `PROVENANCE.txt` are
committed. Every number in §1, §2, and the standardized table in §3 is a cell of one of those three
committed tables. The record-identity check in the provenance section and the fixed-threshold and
metadata-correlation figures in §3 and §4 are derived from `merged.jsonl` directly, by this:

```python
from pathlib import Path
import numpy as np
from motco.simulations.grid import read_replicate_results

recs = read_replicate_results(Path("results/angle-pivotality-2026-09-01/merged.jsonl"))
orient = [r for r in recs
          if (r.cell_metadata or {}).get("trajectory_mode") == "orientation"
          and r.phase == "power_primary"]
null_cell = [r for r in recs
             if r.phase == "type_i_baseline"
             and (r.cell_metadata or {}).get("trajectory_mode") in (None, "none")]

# §3 — one fixed critical value calibrated on the null cell's observed angles
threshold = np.quantile([r.pair_statistics["angle"] for r in null_cell], 0.95)   # 153.22 deg
obs = np.array([r.pair_statistics["angle"] for r in orient])
print(threshold, (obs > threshold).mean())                                       # 0.01 vs 0.68 as specified

# §4 — the null width against the persisted integration diagnostics
q95 = np.log([r.null_summary["angle"]["q95"] for r in orient])
lv = [float(r.integration_metadata["selected_lv"]) for r in orient]
auroc = [float(r.integration_metadata["cv_mean_auroc"]) for r in orient]
print(np.corrcoef(q95, lv)[0, 1], np.corrcoef(q95, auroc)[0, 1])                 # +0.139, -0.179
```
