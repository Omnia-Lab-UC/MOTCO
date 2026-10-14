# Net-displacement Type I study — 2026-10-13

> Paper-grade, null-only re-measurement of the net-displacement orientation estimator
> (`openspec/changes/archive/2026-10-13-net-displacement-type-i-study/`). It settles the Type I question that the
> [orientation measurement probe](orientation-measurement-probe-2026-10-12.md) left open. It runs no power cell, makes no
> power claim, and changes no production default. Every rate below is read from
> `results/net-displacement-type-i-2026-10-13/{phase6,phase5}/report/` (`anchor_decision.csv`/`.json`,
> `design_point_operating.csv`, `type_i_table.csv`, `acceptance_report.csv`) and
> `results/net-displacement-type-i-2026-10-13/follow_up_decision.json`. Rates are rejection rates at p < 0.05 over 500
> replicates × 999 permutations (Monte Carlo SE about 0.0097 at 0.05).

**Runs:** `results/net-displacement-type-i-2026-10-13/phase6/` (SLURM array 917290) and `…/phase5/` (array 917291) on
`ing` · **Code:** `9089729` · **Anchor decision:** `size_holds` at every column at both designs ·
**Follow-up:** `paper_grade_power_comparison`.

**In one paragraph.** The probe found the net-displacement anchor `angle` rate at 0.08–0.10 in every space at both
designs, against 0.01–0.05 under PC1, over 100 replicates. At 500 replicates × 999 permutations that excess does not
reproduce. Every net-displacement column holds the predeclared bound 0.05 + 2·√(0.05·0.95/500) = 0.0695 on `delta`,
`angle` and `shape`, at both designs. The lead arm, `pls_cv × net_displacement`, has anchor `angle` 0.050 at Phase 6
and 0.050 at Phase 5. Net displacement does move the anchor `angle` rate from conservative toward nominal at Phase 6
(PC1 0.024, net displacement 0.048–0.058); it is a size change within the bound, not an excess. The predeclared rule
therefore selects the paper-grade power comparison of `pls_cv × net_displacement` against `pls_cv × pc1`. One caveat
carries forward: under net displacement the translation control's `angle` rate sits near or above the anchor bound in
some columns (Phase 6 `concat` 0.078), so the comparison should hold that control to the same bound (§5).

## 1. Design

- **Columns:** five explicit design points (`design_grid.points`) on matched datasets:

  | column | `integration_method` | `forced_components` | `n_components` | `orientation_estimator` |
  |---|---|---|---|---|
  | `pls_cv × pc1` (baseline) | `pls` | — (CV) | — | `pc1` |
  | `pls_cv × nd` (lead arm) | `pls` | — (CV) | — | `net_displacement` |
  | `pls_forced × nd` | `pls` | 2(K − 1): 4 / 6 | — | `net_displacement` |
  | `pca × nd` | `pca` | — | 2(K − 1): 4 / 6 | `net_displacement` |
  | `concat × nd` (ceiling reference) | `concat` | — | — | `net_displacement` |

  The fixed-basis `pls_within_group` arm is excluded (the probe identified why its anchor inflates).
- **Designs:**
  - **Phase 6** (`examples/trajectory_power_study/net_displacement_type_i_phase6.json`): three stages,
    `group_stage_sizes = [[11, 10, 28], [9, 10, 12]]`, methylation + expression, base seed 800, family `phase6-small-n`.
  - **Phase 5** (`net_displacement_type_i_phase5.json`): n = 1200, four stages, `group_ratio` 0.5, all three blocks, base
    seed 600, family `phase5-primary`.
  - Both use `p_dmp = 0.1`, ρ = 0, joint magnitude and `surgery_censoring = error`.
- **Cells (12 per design):** per column, the shared zero-effect anchor (`none`, e = 0) and the `translation` negative
  control at e = 1.0, plus the engine's two `type_i_baseline` cells (independently seeded, measured in the baseline
  column).
- **Seeds:** every column evaluates the same datasets at each replicate index. `evaluation.seed` is the constant 3 and is
  also the RRPP seed, so every column draws the same permutation indices, and the baseline column re-computes the parent
  studies' records (§4).
- **Rule:** `acceptance.anchor_decision` (rule version 1): per column and statistic, anchor rate ≤ α + 2·√(α(1 − α)/R)
  with α = 0.05 and R the available replicates (0.0695 at R = 500). A column is `size_holds` when all three statistics
  hold. Advisory; it feeds neither the gate (off) nor the acceptance targets.

## 2. Anchor rates against the bound

From `anchor_decision.csv`. Rate (MC SE); bound 0.0695 in every cell (R = 500 everywhere, no failed replicate).

**Phase 6 (n = 80, K = 3).**

| column | `delta` | `angle` | `shape` | verdict |
|---|---|---|---|---|
| `pls_cv × pc1` | 0.026 (0.0071) | 0.024 (0.0068) | 0.020 (0.0063) | `size_holds` |
| `pls_cv × nd` | 0.026 (0.0071) | **0.050** (0.0097) | 0.020 (0.0063) | `size_holds` |
| `pls_forced4 × nd` | 0.006 (0.0035) | 0.048 (0.0096) | 0.006 (0.0035) | `size_holds` |
| `pca4 × nd` | 0.004 (0.0028) | 0.054 (0.0101) | 0.000 (0.0000) | `size_holds` |
| `concat × nd` | 0.006 (0.0035) | 0.058 (0.0105) | 0.002 (0.0020) | `size_holds` |

**Phase 5 (n = 1200, K = 4).**

| column | `delta` | `angle` | `shape` | verdict |
|---|---|---|---|---|
| `pls_cv × pc1` | 0.036 (0.0083) | 0.046 (0.0094) | 0.048 (0.0096) | `size_holds` |
| `pls_cv × nd` | 0.036 (0.0083) | **0.050** (0.0097) | 0.048 (0.0096) | `size_holds` |
| `pls_forced6 × nd` | 0.016 (0.0056) | 0.048 (0.0096) | 0.000 (0.0000) | `size_holds` |
| `pca6 × nd` | 0.008 (0.0040) | 0.052 (0.0099) | 0.006 (0.0035) | `size_holds` |
| `concat × nd` | 0.008 (0.0040) | 0.050 (0.0097) | 0.006 (0.0035) | `size_holds` |

`delta` and `shape` do not depend on the estimator: in the two `pls_cv` columns they are the same records' statistics
and their rates are identical, as they must be. The fixed-rank and label-free spaces are conservative on `delta` and
`shape` at both designs; that is a property of those spaces, not of the estimator, and no rule here reads it.

**Against the probe.** The probe's lead-arm anchor `angle` rates were 0.100 (Phase 6) and 0.080 (Phase 5) over
100 replicates × 199 permutations, on its own seed root (base seed 900), so its datasets are independent of these.
Here the rate is 0.050 at both designs. Split into five consecutive blocks of 100 replicates, the lead arm's Phase 6
anchor `angle` rate ranges from 0.01 to 0.09 (0.06, 0.09, 0.01, 0.04, 0.05), and at Phase 5 from 0.02 to 0.07 (0.07,
0.04, 0.02, 0.05, 0.07); a 100-replicate rate of 0.08–0.10 is within that sampling spread. The probe's excess was
Monte Carlo noise at R = 100, which is what the probe report said it could not exclude.

## 3. Negative controls (descriptive)

From `anchor_decision.csv` (`role = control`). Translation at e = 1.0 moves group B's location only, so it is a second
null for `angle` and `shape`; the rule does not classify it.

| design | column | translation `delta` | translation `angle` | translation `shape` |
|---|---|---|---|---|
| Phase 6 | `pls_cv × pc1` | 0.028 | 0.028 | 0.030 |
| Phase 6 | `pls_cv × nd` | 0.028 | 0.068 | 0.030 |
| Phase 6 | `pls_forced4 × nd` | 0.002 | 0.058 | 0.000 |
| Phase 6 | `pca4 × nd` | 0.004 | 0.064 | 0.000 |
| Phase 6 | `concat × nd` | 0.004 | **0.078** | 0.002 |
| Phase 5 | `pls_cv × pc1` | 0.042 | 0.042 | 0.030 |
| Phase 5 | `pls_cv × nd` | 0.042 | 0.040 | 0.030 |
| Phase 5 | `pls_forced6 × nd` | 0.006 | 0.064 | 0.004 |
| Phase 5 | `pca6 × nd` | 0.012 | 0.064 | 0.000 |
| Phase 5 | `concat × nd` | 0.012 | 0.064 | 0.000 |

The `type_i_baseline` cells (baseline column, independent datasets) give, at Phase 6, `none` 0.016 / 0.028 / 0.018
and `translation` 0.034 / 0.040 / 0.018 (`delta` / `angle` / `shape`); at Phase 5, `none` 0.032 / 0.044 / 0.050 and
`translation` 0.036 / 0.056 / 0.070 (`type_i_table.csv`). The two-sided acceptance targets flag the conservative
Phase 6 rates (and Phase 5 `none`/`delta` 0.032) as failures; they are below α, not above it, which is the two-sided
target's known behavior (roadmap, next-change item 2). Phase 6's two Type I cells are the parent study's cells and
reproduce its records exactly (§4), so those rates are already on record there.

## 4. Covariates and reproduction of the parent studies

**Eigengap and `angle` null width** (`design_point_operating.csv`, anchor rows; median recorded pooled eigengap, median
and IQR of the per-replicate `angle` null q95, in degrees):

| design | column | median eigengap | `angle` null q95 median (IQR) | selected rank |
|---|---|---|---|---|
| Phase 6 | `pls_cv × pc1` | 0.111 | 84.0° (115.4) | CV, median 2 (2–3) |
| Phase 6 | `pls_cv × nd` | 0.111 | 12.1° (26.6) | CV, median 2 (2–3) |
| Phase 6 | `pls_forced4 × nd` | 0.133 | 32.8° (5.9) | forced 4 |
| Phase 6 | `pca4 × nd` | 0.107 | 38.3° (3.6) | — |
| Phase 6 | `concat × nd` | 0.105 | 43.8° (3.4) | — |
| Phase 5 | `pls_cv × pc1` | 0.053 | 5.5° (4.9) | CV, median 3 (2–3) |
| Phase 5 | `pls_cv × nd` | 0.053 | 0.7° (0.3) | CV, median 3 (2–3) |
| Phase 5 | `pls_forced6 × nd` | 0.051 | 9.2° (0.9) | forced 6 |
| Phase 5 | `pca6 × nd` | 0.052 | 9.5° (0.9) | — |
| Phase 5 | `concat × nd` | 0.052 | 11.8° (1.1) | — |

Net displacement narrows the `angle` null by a factor of 7–8 in the `pls_cv` space (Phase 6 84.0° → 12.1°, Phase 5 5.5° → 0.7°) and removes
most of its across-replicate dispersion (Phase 6 IQR 115° → 27°), which is the PC1 sign/flip tail the probe described.
A narrower null is consistent with size holding: the test calibrates against its own permutation distribution.

**Reproduction.** The baseline column shares the parent study's generator seeds and evaluation parameters, including
the RRPP seed, so it re-computes the parent's records:

- **Phase 6:** all 1,000 baseline records (anchor and translation, replicates 0–499) equal
  `results/phase6-small-n-2026-10-09/merged.jsonl` in generator seed, p-values and observed statistics. The two
  `type_i_baseline` cells carry the parent's cell ids and reproduce its 1,000 records' p-values exactly.
- **Phase 5:** all 1,000 baseline records equal `results/phase5-2026-09-10/merged.jsonl` in generator seed, selected rank
  and p-values. Observed statistics are bit-identical in 317 records and differ in the rest by at most 2.8 × 10⁻¹⁰
  relative (floating point; the parent ran at `b8f7daa`), never enough to change a p-value. The `type_i_baseline`
  cells are new (`magnitude_kind` enters their ids).

The baseline anchor rates therefore equal the parent studies' recorded rates (Phase 6 0.026 / 0.024 / 0.020; Phase 5
0.036 / 0.046 / 0.048), and the five columns are paired measurements of the same 500 datasets per design.

## 5. Decision and the next change

`anchor_decision.json`, both designs: `size_holds` at all five columns. `follow_up_decision.json` (rule version 1, inputs
`phase6/report/anchor_decision.json` sha256 `a0af9d20…`, `phase5/report/anchor_decision.json` sha256 `89e8e449…`):

- **verdict `paper_grade_power_comparison`**: the lead arm `pls_cv × net_displacement` is `size_holds` at both
  designs;
- localization `none` at both designs (no net-displacement column exceeds).

The verdict selects the next change: **a paper-grade power comparison of `pls_cv × net_displacement` against
`pls_cv × pc1`** at both designs, through the study engine (the estimator is now a contracted design axis, so the
comparison is one explicit-points config per design over the Phase 5 / Phase 6 power grids). The anchor-excess
diagnosis is not needed. Two points to carry into that change:

1. **Hold the translation control to the anchor bound.** Under net displacement the translation `angle` rate is
   0.068 for the lead arm at Phase 6 (bound 0.0695) and 0.058–0.078 in the other spaces; Phase 6 `concat` (0.078)
   would exceed the bound if it were judged by it. The lead arm passes at both designs, but by one rejection at
   Phase 6, so the comparison should declare translation `angle` as a mandatory control rather than leave it
   descriptive.
2. **`concat` stays the ceiling reference, not a candidate**, and its p-padding question remains open; this study's
   null-only design cannot answer it.

## 6. Provenance and reproduction

- **Configs:** `examples/trajectory_power_study/net_displacement_type_i_phase6.json` (sha256 `559ca569…`) and
  `…_phase5.json` (sha256 `304838d8…`), identity-pinned in `tests/data/committed_config_identity.json`.
- **Code:** `9089729fb071eb291ac0de0300a94ceb0b473839` (branch `feat/net-displacement-type-i-study`), clean tree on the
  cluster. Python 3.11.16, motco 0.6.0, numpy 2.3.5, scikit-learn 1.8.0, scipy 1.16.3, pandas 2.3.3, uv 0.12.10.
- **Host:** `ing`, partition `512x1024` (AMD EPYC 7662), 50 single-CPU shards per design, BLAS pinned to one thread,
  `STUDY_N_JOBS` unset, nodes n4 (memory) and n11 (NFS stalls in the timing job) excluded.
- **Records:** 6,000/6,000 units per design, 0 failed, 0 censored, every signature equal to the enumerated cell's.
- **Cost:** Phase 6 35.2 recorded core-hours (1 h 07 m wall), Phase 5 97.9 (2 h 43 m wall), against timing-job
  projections of 31 and 97 with margin (`PROVENANCE.txt`).
- **Merged records:** `phase6/merged.jsonl` sha256 `fef16757…`, `phase5/merged.jsonl` sha256 `641020f0…` (gitignored;
  verified identical on the cluster and the workstation).

Reproduce from the shards:

```bash
B=results/net-displacement-type-i-2026-10-13
for d in phase6 phase5; do
  python scripts/motco_study.py merge --out-dir $B/$d
  python scripts/motco_study.py report \
    --config examples/trajectory_power_study/net_displacement_type_i_$d.json --out-dir $B/$d
done
python scripts/net_displacement_follow_up.py \
  --phase6 $B/phase6/report/anchor_decision.json \
  --phase5 $B/phase5/report/anchor_decision.json \
  --out $B/follow_up_decision.json
```

The SLURM launch lines, the timing job, and unit timings are in `phase6/PROVENANCE.txt` and `phase5/PROVENANCE.txt`.
