# The recorded eigengap reproduces the audit's association — and stratifies the orientation shortfall

**Run date:** 2026-09-02
**Code revision:** working tree of the `record-latent-config-spectrum` change, parent
`84061bdec824da9c33865d52e59a13b3750fb9fb` (the merge of the geometry audit)
**Configuration:** `examples/trajectory_power_study/angle_pivotality_diagnostic.json` — byte-identical to
the 2026-09-01 diagnostic run (same SHA-256, `2dfa5cc3…`)
**Records:** `results/latent-spectrum-2026-09-02/merged.jsonl` (500 replicates, 0 failures, 199 permutations
each, 500 of 500 carrying the configuration-spectrum block)
**Tables:** `results/latent-spectrum-2026-09-02/report/` — `pivotality_spectrum.csv`,
`eigengap_stratified_power.csv`, `config_spectrum.csv`, plus the pre-existing pivotality and study tables
**Provenance:** `results/latent-spectrum-2026-09-02/PROVENANCE.txt`
**Change:** `openspec/changes/record-latent-config-spectrum`
**Answers:** [Phase 5 readiness](https://github.com/Omnia-Lab-UC/MOTCO/blob/c440070dcb3c5bcc6ddd2227cf7adbef18573b6a/docs/phase5-readiness.md) items 1 (left-open gap), 3, and 4 (stratifying
covariate); geometry audit plan item **P1**
**Environment:** Python 3.11.15, numpy 2.3.5, pandas 2.3.3, scikit-learn 1.8.0, scipy 1.16.3; no R at runtime

## Summary

The [pivotality report](angle-null-pivotality-2026-09-01.md) established that an orientation replicate's
outcome is decided by the width of its own `angle` permutation null, and stated that *"nothing the harness
records today identifies which replicates will be resolvable"*. The [geometry
audit](geometry-audit-2026-09-01.md) (finding F1) named the missing quantity — the relative eigengap
`(λ₁−λ₂)/Σλ` of the centered latent stage-mean configuration — by regenerating replicates from persisted
seeds. This run **persists that quantity per replicate** and re-measures the association from the record
alone.

Three results, all from the same 500 records:

1. **The audit's association reproduces.** In the orientation cell, Spearman(pooled eigengap, own `angle`
   null q95) = **−0.749** [−0.824, −0.648]; log–log Pearson **−0.622**. The audit measured −0.75 and −0.62.
2. **It is specific to `angle`.** In the same cell the association with the `delta` and `shape` null widths
   is +0.093 [−0.105, +0.284] and −0.032 [−0.227, +0.166] — both intervals cover zero.
3. **It stratifies the shortfall.** Orientation `angle` power rises monotonically across within-cell
   eigengap terciles: **0.44 → 0.76 → 0.85**. The overall 0.68 is an average over a covariate the harness
   now records.

And the change is inert: every rejection rate and every pivotality slope reproduces the 2026-09-01 run
exactly.

## 1. The association, measured from records

`pivotality_spectrum.csv`, `null_target = q95`, n = 100 per cell. Spearman with a Fisher-z 95% interval:

| Cell | `delta` | `angle` | `shape` |
|---|---|---|---|
| `none` (zero effect) | −0.092 | **−0.814** [−0.871, −0.735] | −0.108 |
| `translation` @1.00 (`type_i_baseline`) | +0.019 | **−0.740** [−0.817, −0.636] | +0.163 |
| `translation` @1.00 (`power_primary`) | +0.053 | **−0.836** [−0.887, −0.765] | +0.009 |
| `magnitude` @1.00 | −0.086 | **−0.701** [−0.789, −0.586] | −0.038 |
| `orientation` @1.00 | +0.093 | **−0.749** [−0.824, −0.648] | −0.032 |

Every `angle` interval excludes zero by a wide margin; no `delta` or `shape` interval does. Log–log Pearson
tells the same story for `angle` (−0.61 to −0.71 across cells; −0.622 in the orientation cell).

**Read this as a property of the geometry, not of the effect.** The association is present in the zero-effect
baseline and in both translation controls at essentially the same strength. That is what the audit's
mechanism predicts: PC1's estimator variance scales like noise over the eigengap, so a near-isotropic stage
configuration yields a wide `angle` null whether or not a group difference was constructed. The recorded
eigengap therefore predicts *resolvability*, and it does so without seeing the group labels.

## 2. It stratifies the orientation shortfall

`eigengap_stratified_power.csv`, orientation @1.00, α = 0.05, within-cell terciles of the recorded pooled
eigengap:

| Stratum | eigengap range | mean eigengap | n | `angle` rejections | rate | MC SE |
|---|---|---|---|---|---|---|
| bottom | 0.0057 – 0.0359 | 0.0243 | 34 | 15 | **0.44** | 0.085 |
| middle | 0.0364 – 0.0552 | 0.0451 | 33 | 25 | **0.76** | 0.075 |
| top | 0.0556 – 0.1173 | 0.0739 | 33 | 28 | **0.85** | 0.062 |

A 3× separation in eigengap maps to a 0.41 separation in power (SE of the difference 0.105). The cell's
headline 0.68 is an average over this spread: replicates whose latent stage configuration has a dominant
axis clear the 0.80 floor; those whose configuration is near-isotropic do not, and no amount of signal in a
direction the configuration cannot resolve will move them.

The recorded distributions (`config_spectrum.csv`) show how tight the regime is. Pooled eigengaps per cell:

| Cell | mean | median | min | max |
|---|---|---|---|---|
| `orientation` @1.00 | 0.0475 | 0.0444 | 0.0057 | 0.1173 |
| `magnitude` @1.00 | 0.0422 | 0.0406 | 0.0028 | 0.1612 |
| `translation` @1.00 (`power_primary`) | 0.0431 | 0.0404 | 0.0080 | 0.1126 |
| `translation` @1.00 (`type_i_baseline`) | 0.0532 | 0.0537 | 0.0030 | 0.1670 |
| `none` | 0.0450 | 0.0409 | 0.0022 | 0.1248 |

Every cell sits near the isotropic end (a straight trajectory would read 1.0), and the per-group entries
agree with the pooled ones — consistent with the audit's account that the independent-indicator baseline,
not the trajectory surgery, sets the eigengap. That is the lever readiness item 4 hands to baseline
continuity (plan item P4), and it is now measured per replicate rather than inferred.

## 3. The change is inert

The config is byte-identical to the 2026-09-01 diagnostic, so the operating point is directly comparable.
Every published number reproduces:

| Quantity | 2026-09-01 | this run |
|---|---|---|
| orientation `angle` / `delta` / `shape` power | 0.68 / 0.65 / 0.99 | **0.68 / 0.65 / 0.99** |
| magnitude `delta` / `angle` / `shape` | 1.00 / 0.00 / 0.00 | **1.00 / 0.00 / 0.00** |
| Type I `none` (delta/angle/shape, combined) | — | 0.04 / 0.05 / 0.01, 0.09 |
| `angle` null-tracking slope: `none` / `translation` (t1b, pp) / `magnitude` / `orientation` | 0.931 / 0.925, 0.960 / 0.874 / 0.811 | **0.931 / 0.925, 0.960 / 0.874 / 0.811** |
| orientation `delta` / `shape` slope | 0.018 / 0.058 | **0.018 / 0.058** |

Reproduction to three decimals on every slope, and to the replicate on every rate. The spectrum consumes no
randomness: it is computed from the LS-mean vectors `estimate_difference` already fits, after the
permutation index is drawn. Unit tests pin the permutation draws and null distributions as **byte-identical**
with the recording flag on and off, serially and in parallel
(`tests/test_trajectory_spectrum.py`).

**Cost.** One SVD of a k×d matrix per configuration. On this design point (4 stages, rank-3 PLS latent
space) recording the per-permutation eigengap adds ~12 ms to a 199-permutation RRPP — under 0.1% of a ~14 s
replicate, which is dominated by the PLS double cross-validation.

## 4. What this does and does not change

- **Nothing statistical.** No statistic, p-value, permutation draw, generator behaviour, integration method,
  or latent-space sizing changed. The eigengap is a recorded covariate and a reporting qualifier — never a
  filter, weight, or test modification.
- **The audit's caution stands.** The spectrum is only *observed*. Any group-aware sizing or supervision
  would void the fixed-latent-space RRPP conditioning; nothing here reads the group labels.
- **Resume policy changed deliberately.** `parameter_signature` now includes `config_spectrum_version`, so a
  shard written before this change refuses resume rather than producing a merged set in which only some
  records carry the covariate. Already-merged pre-change record sets still load (empty block) and their
  committed reports are unaffected — verified by reporting a pre-change record set under the new code and
  finding `specificity_matrix.csv`, `power_curves.csv`, `type_i_table.csv`, `acceptance_report.csv`, and
  `acceptance_report.json` byte-identical to the pre-change code's output.

## Reproduction

```bash
# 5 cells x 100 replicates = 500 work units; ~50 minutes on 16 local shards.
for i in $(seq 0 15); do
  uv run python scripts/run_study_shard.py \
    --config examples/trajectory_power_study/angle_pivotality_diagnostic.json \
    --out-dir results/latent-spectrum-2026-09-02 \
    --shard-index "$i" --n-shards 16 --error-policy record &
done
wait
uv run python scripts/motco_study.py merge --out-dir results/latent-spectrum-2026-09-02
uv run python scripts/motco_study.py report \
  --config examples/trajectory_power_study/angle_pivotality_diagnostic.json \
  --out-dir results/latent-spectrum-2026-09-02
uv run python scripts/angle_null_pivotality.py \
  --merged  results/latent-spectrum-2026-09-02/merged.jsonl \
  --out-dir results/latent-spectrum-2026-09-02/report
```

Do **not** pass `--n-jobs`: it is part of the cell parameter signature. The shard and merged JSONL are
gitignored as regenerable; `report/` and `PROVENANCE.txt` are committed.
