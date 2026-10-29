# Design-guidance paper-grade slice — findings report (2026-10-29)

**Change:** `design-guidance-paper-grade-slice` (MOTCO 0.12.0). **Supersedes at its columns:** the pilot-grade rates of the [design-guidance pilot](design-guidance-pilot-2026-10-28.md) (2026-10-28) at n ∈ {80, 200, 500}, three stages, two or three blocks, and at the joint stage effect δ ∈ {0.25, 0.5, 1.0}; the pilot's four-stage, n = 1000 and remaining arm columns stand at pilot grade. **Reads:** `results/design-guidance-slice-2026-10-29/` (`report/` and `PROVENANCE.txt` per profile, `figures/`, `reproduction/`).

The slice reruns the pilot report §7 proposal at paper grade: the three-stage core at the sample sizes the planned analyses have, under the production measurement (pooled PLS with the stage-supervised double-CV rank floored at 2, net-displacement orientation, joint magnitude), 500 replicates × 999 permutations (SE ≤ 0.022 on every rate; anchor bound 0.0695 at R = 500), in the pilot's seed family so that every column evaluates the pilot's datasets at its first 100 replicates. It is guidance, not a gate: the predeclared rules are the advisory anchor decision and a column-scoped Type I rule for translation's `angle`; there is no power target. The second-priority stage-effect arm (§8) was run on the user's decision once the slice's numbers were in.

**In one paragraph.** Size holds at every column under both predeclared rules. `delta` is the component a small cohort detects: 0.62 at n = 80 and 0.94 at n = 200 at e = 0.1 with two blocks (0.84 and 0.99 with three), and 0.63 at n = 500 even at e = 0.02 with two blocks. `shape` is monotone in n, effect and block count: ≤ 0.20 at n = 80 with two blocks, 0.56 at n = 200 and 0.85 at n = 500 at e = 1.0, and 0.48 / 0.86 / 0.97 with three blocks. `angle` is governed by n and block count with a flat, non-monotone response in effect (0.10–0.15 at n = 80 and 0.16–0.24 at n = 200 with two blocks; 0.36–0.43 at n = 500; 0.24–0.31 / 0.41–0.49 / 0.66–0.74 with three), every reading at a CV-selected rank of 2 or 3. All of it at a stage signal that is saturated (held-out CV AUROC ≥ 0.9998 at every cell), which is what the stage-effect arm relaxes.

## 1. Configuration and provenance

| | core magnitude | core geometry | stage-effect magnitude | stage-effect geometry |
|---|---|---|---|---|
| config | `design_guidance_core_magnitude_study.json` | `design_guidance_core_geometry_study.json` | `design_guidance_stage_effect_magnitude_study.json` | `design_guidance_stage_effect_geometry_study.json` |
| config sha256 (first 16) | `3f4a618d57fc59db` | `bb56f1c6f61ce2e3` | `2cd406d7d6238ff0` | `961bac86a2e2a6c8` |
| design columns | n {80, 200, 500} × 3 stages × {2, 3} blocks (6) | same (6) | joint δ {0.25, 0.5, 1.0} at n = 200, 3 stages, 2 blocks (3) | same (3) |
| modes × effects | `magnitude` × {0.02, 0.05, 0.1} | `orientation`, `shape`, `translation` × {0.25, 0.5, 1.0} | `magnitude` × {0.02, 0.05, 0.1} | `orientation`, `shape`, `translation` × {0.25, 0.5, 1.0} |
| cells / units | 26 / 13,000 | 62 / 31,000 | 14 / 7,000 | 32 / 16,000 |
| SLURM jobs (shards) | 929844 + 929980 (40) | 929845 (95) | 930041 (20) | 930042 (46) |

All four profiles: 500 replicates × 999 permutations, `base_seed` 1028, matched-seed family `design-guidance`, `evaluation.seed` 3 (also the RRPP seed), `orientation_estimator` `net_displacement`, `min_components` 2, `magnitude_kind` `joint`, the Phase 6 CV knobs (cv1 3, cv2 4, 5 repeats, ≤ 20 components), attribution and gate off, `surgery_censoring` `error`. The slice baseline is the pilot's (n = 200, three stages, `group_ratio` 0.5, `p_dmp` 0.1, ρ = 0, every per-omic δ 2.0, methylation + expression measured); the arm's configuration baseline is δ = 0.25 (§8). Code revision **`8bf8eb7`** on `feat/design-guidance-paper-grade-slice` (clean worktree `/home1/tgonza/MOTCO-design-guidance-slice` on `ing`, the main checkout's `.venv` with `PYTHONPATH` on the worktree's `src`; the venv's `motco` distribution metadata says 0.11.0 and the revision is the version of record); python 3.11.16, numpy 2.3.5, scikit-learn 1.8.0, scipy 1.16.3. Partition `512x1024`, 1 CPU and 4 GB per task, 6 h limit, BLAS pinned to one thread, `STUDY_N_JOBS` unset, `--error-policy record`. Every record carries `n_jobs` 1.

Pre-flight: the fast suite at the same tree (ruff, mypy, 1583 passed / 1 skipped); the change's six test modules on `ing` from the worktree (115 passed); a 2 × 19 smoke of all four profiles (268 / 268 units, 0 failures); 14 units timed at 999 permutations on `ing`, on a quiet node (8–19 s per unit across the six columns) and on a saturated one (16–33 s) — 999 permutations cost no more than the pilot's 199 did, because the PLS double CV dominates a unit. Shards were sized at the saturated cost plus 20 % at ≤ 2.5 h each. Thirty-five magnitude shards that SLURM placed on an overloaded node sat in NFS wait for 32 minutes without writing a unit; they were cancelled and the same array ids resubmitted elsewhere (job 929980), which changes nothing in the records (shards are keyed by the `(cell_id, replicate_index)` hash and are resumable). Full details: `openspec/changes/design-guidance-paper-grade-slice/notes.md`.

## 2. Unit, failure and reproduction accounting

| Profile | Array tasks | Records (expected) | Duplicates / missing / signature mismatches | Failed | Censored | Recorded core-hours (median unit) | Wall |
|---|---|---:|---|---:|---:|---:|---|
| core magnitude | 40 `COMPLETED 0:0` (+ 35 cancelled and resubmitted before writing) | 13,000 (13,000) | 0 / 0 / 0 | 0 | 0 | 52.7 (11.9 s) | 11:42 → 14:27 |
| core geometry | 95 `COMPLETED 0:0` | 31,000 (31,000) | 0 / 0 / 0 | 0 | 0 | 119.0 (12.1 s) | 11:42 → 13:29 |

Every record carries `n_jobs` 1, `orientation_estimator` `net_displacement` and `min_components` 2. The slice's 44,000 units cost 171.7 recorded core-hours, between the quiet-node (151) and saturated-node (278) estimates and below the design's 380. The CV selected rank 2 in 35,000 of the 44,000 units and rank 3 in the other 9,000; never more.

**Reproduction** (`scripts/anchor_reproduction.py --by-column`, `reproduction/`):

- *Siblings* (`siblings.csv`): the magnitude and geometry profiles' zero-effect anchors pair at every one of the six columns, 500 of 500 replicates each, 3,000 pairs: generator seeds, selected ranks, observed `delta`/`angle`/`shape` (max |Δ| ≤ 5 × 10⁻¹¹) and **p-values identical** on all 3,000. The report reads the geometry profile's anchor; the magnitude profile's is this check.
- *Pilot* (`pilot_core_magnitude.csv`, `pilot_core_geometry.csv`): against the pilot core profiles' `merged.jsonl` (199 permutations), every slice column pairs replicate indices 0–99 — 600 pairs per profile — with seeds, ranks and observed statistics identical (max |Δ| ≤ 4 × 10⁻¹⁰); p-values are reported as not comparable and excluded from the verdict; the pilot's ten other columns (four stages; n = 1000) are listed with zero pairs. The slice therefore evaluates the pilot's datasets at its first 100 replicates and 400 fresh ones; it supersedes the pilot's rates at these columns and never pools the two.

## 3. Per-column operating characteristics

Rates are rejection rates at α = 0.05 over 500 available replicates, ± 1 Monte Carlo SE. The figures are `figures/design_guidance_{magnitude,orientation,shape,control}.png` (`scripts/design_guidance_figure.py --grade paper`, one panel per block count at three stages, n on a log axis, one line per effect) with the plotted rows in `figures/design_guidance_power.csv`; the full tables, with the stage signal and selected rank per cell, are `results/design-guidance-slice-2026-10-29/tables.md`.

### 3.1 Magnitude / `delta`

| column | e = 0.02 | e = 0.05 | e = 0.1 |
|---|---|---|---|
| n = 80, 2 blocks | 0.036 ± 0.008 | 0.234 ± 0.019 | 0.620 ± 0.022 |
| n = 80, 3 blocks | 0.112 ± 0.014 | 0.512 ± 0.022 | 0.836 ± 0.017 |
| n = 200, 2 blocks | 0.228 ± 0.019 | 0.694 ± 0.021 | 0.936 ± 0.011 |
| n = 200, 3 blocks | 0.430 ± 0.022 | 0.878 ± 0.015 | 0.992 ± 0.004 |
| n = 500, 2 blocks | 0.632 ± 0.022 | 0.962 ± 0.009 | 0.998 ± 0.002 |
| n = 500, 3 blocks | 0.834 ± 0.017 | 0.994 ± 0.003 | 1.000 ± 0.000 |

Monotone in n, effect and block count at every column. The pilot's values at these 18 cells are reproduced within 0.04 (mean difference +0.004), inside the pilot's SE of ≤ 0.05.

### 3.2 Orientation / `angle` and shape / `shape`

| column | `angle` e = 0.25 | `angle` e = 0.5 | `angle` e = 1.0 | `shape` e = 0.25 | `shape` e = 0.5 | `shape` e = 1.0 |
|---|---|---|---|---|---|---|
| n = 80, 2 blocks | 0.100 ± 0.013 (rank 3) | 0.110 ± 0.014 (3) | 0.150 ± 0.016 (3) | 0.044 ± 0.009 | 0.104 ± 0.014 | 0.204 ± 0.018 |
| n = 80, 3 blocks | 0.238 ± 0.019 (2) | 0.312 ± 0.021 (2) | 0.274 ± 0.020 (3) | 0.184 ± 0.017 | 0.316 ± 0.021 | 0.484 ± 0.022 |
| n = 200, 2 blocks | 0.158 ± 0.016 (2) | 0.242 ± 0.019 (3) | 0.184 ± 0.017 (3) | 0.236 ± 0.019 | 0.418 ± 0.022 | 0.562 ± 0.022 |
| n = 200, 3 blocks | 0.406 ± 0.022 (2) | 0.494 ± 0.022 (2) | 0.492 ± 0.022 (2) | 0.486 ± 0.022 | 0.668 ± 0.021 | 0.856 ± 0.016 |
| n = 500, 2 blocks | 0.364 ± 0.022 (2) | 0.426 ± 0.022 (2) | 0.358 ± 0.021 (3) | 0.548 ± 0.022 | 0.750 ± 0.019 | 0.848 ± 0.016 |
| n = 500, 3 blocks | 0.662 ± 0.021 (2) | 0.738 ± 0.020 (2) | 0.710 ± 0.020 (2) | 0.780 ± 0.019 | 0.892 ± 0.014 | 0.970 ± 0.008 |

The number in parentheses is the median CV-selected rank at that cell (the range is 2–3 at every cell of the slice). `shape` is monotone in n, effect and block count. `angle` rises with n and with the third block but not with the effect: at every column the e = 1.0 rate is within 0.08 of the e = 0.25 rate and below the e = 0.5 rate at four of six columns, the pattern the pilot saw at pilot grade, now resolved at SE ≤ 0.022. Against the pilot the orientation rates are reproduced within 0.11 (mean +0.044, the largest differences at n = 200 and 500 with three blocks at e = 1.0, where the pilot's 100 replicates read 0.38 and 0.61 and the slice reads 0.49 and 0.71) and the shape rates within 0.12 (mean +0.016).

#### Translation (negative control), all statistics at every effect

| column | e | `delta` | `angle` | `shape` | column-scoped verdict (`angle` at e = 1.0 vs 0.0695) |
|---|---|---|---|---|---|
| n = 80, 2 blocks | 0.25 | 0.020 ± 0.006 | 0.074 ± 0.012 | 0.020 ± 0.006 | |
| | 0.5 | 0.034 ± 0.008 | 0.078 ± 0.012 | 0.020 ± 0.006 | |
| | 1.0 | 0.032 ± 0.008 | 0.040 ± 0.009 | 0.024 ± 0.007 | **holds** |
| n = 80, 3 blocks | 0.25 | 0.044 ± 0.009 | 0.064 ± 0.011 | 0.062 ± 0.011 | |
| | 0.5 | 0.058 ± 0.010 | 0.056 ± 0.010 | 0.044 ± 0.009 | |
| | 1.0 | 0.046 ± 0.009 | 0.028 ± 0.007 | 0.048 ± 0.010 | **holds** |
| n = 200, 2 blocks | 0.25 | 0.042 ± 0.009 | 0.054 ± 0.010 | 0.040 ± 0.009 | |
| | 0.5 | 0.058 ± 0.010 | 0.054 ± 0.010 | 0.050 ± 0.010 | |
| | 1.0 | 0.050 ± 0.010 | 0.052 ± 0.010 | 0.032 ± 0.008 | **holds** |
| n = 200, 3 blocks | 0.25 | 0.046 ± 0.009 | 0.052 ± 0.010 | 0.042 ± 0.009 | |
| | 0.5 | 0.062 ± 0.011 | 0.050 ± 0.010 | 0.056 ± 0.010 | |
| | 1.0 | 0.048 ± 0.010 | 0.044 ± 0.009 | 0.050 ± 0.010 | **holds** |
| n = 500, 2 blocks | 0.25 | 0.042 ± 0.009 | 0.054 ± 0.010 | 0.068 ± 0.011 | |
| | 0.5 | 0.042 ± 0.009 | 0.048 ± 0.010 | 0.050 ± 0.010 | |
| | 1.0 | 0.056 ± 0.010 | 0.048 ± 0.010 | 0.042 ± 0.009 | **holds** |
| n = 500, 3 blocks | 0.25 | 0.050 ± 0.010 | 0.048 ± 0.010 | 0.054 ± 0.010 | |
| | 0.5 | 0.036 ± 0.008 | 0.046 ± 0.009 | 0.052 ± 0.010 | |
| | 1.0 | 0.050 ± 0.010 | 0.032 ± 0.008 | 0.052 ± 0.010 | **holds** |

Translation runs at all three geometry effects because a configuration has one effect list (design D2); the predeclared rule reads e = 1.0 and the two smaller offsets are descriptive. No translation rate exceeds 0.078 anywhere. The two readings above the bound are `angle` at the two smaller offsets at n = 80 with two blocks (0.074 and 0.078, SE 0.012); at the same column the e = 1.0 rate is 0.040 and the anchor's 0.036.

## 4. Stage signal and selected rank

The held-out CV stage AUROC (`design_point_stage_signal.csv`, the number that compares with a real analysis) is ≥ 0.9998 in mean and 1.000 at the q33–q67 range at every cell of both profiles: the stage signal is saturated at the core stage effect δ = 2.0, as in the pilot. The CV selects rank 2 or 3 at every cell (median 2 at 36 of the 42 geometry cells and 3 at the others, the n = 80 / two-block column and the e ≥ 0.5 orientation cells at n = 200 and 500 with two blocks among them; never above 3). With a saturated stage signal and two or three retained components, orientation is measured in a two- or three-dimensional latent space whose second axis carries little signal, and `angle` power is governed by how much of the group difference that second axis resolves — by n and the block count — rather than by the size of the orientation surgery (§9).

## 5. Size at every column

**Anchor decision** (`anchor_decision.json`; zero-effect anchor, R = 500, bound 0.0695; the geometry profile's anchor, the magnitude profile's being its bit-exact reproduction):

| column | `delta` | `angle` | `shape` | verdict |
|---|---|---|---|---|
| n = 80, 2 blocks | 0.022 ± 0.007 | 0.036 ± 0.008 | 0.010 ± 0.004 | size_holds |
| n = 80, 3 blocks | 0.052 ± 0.010 | 0.032 ± 0.008 | 0.046 ± 0.009 | size_holds |
| n = 200, 2 blocks | 0.040 ± 0.009 | 0.060 ± 0.011 | 0.034 ± 0.008 | size_holds |
| n = 200, 3 blocks | 0.038 ± 0.009 | 0.060 ± 0.011 | 0.062 ± 0.011 | size_holds |
| n = 500, 2 blocks | 0.040 ± 0.009 | 0.050 ± 0.010 | 0.040 ± 0.009 | size_holds |
| n = 500, 3 blocks | 0.042 ± 0.009 | 0.046 ± 0.009 | 0.040 ± 0.009 | size_holds |

Every statistic holds at every column; the largest anchor rate is 0.062 (`shape`, n = 200, three blocks) against 0.0695. Chance predicts about 0.4 excesses over these 18 tests; none occurred. The pilot's two excesses (both 0.10 against 0.0936 at R = 100) were at columns outside the slice.

**Column-scoped Type I rule** (`acceptance_report.csv`, `type_i_control[translation@1.0, angle, <column>]`, upper bound at α 0.05 and k 2): translation's `angle` at e = 1.0 is **met at all six columns** of the geometry profile, at 0.028–0.052 against 0.0695 (the table in §3.2). In the magnitude profile the same entry reports every column unavailable, as predeclared: translation runs only in its Type I baseline phase, and the entry is declared so that the sibling acceptance blocks are identical.

**Type I baseline phase** (the two cells keyed on their own ids, a fresh 500-replicate draw at the baseline column, two-sided entry at tolerance 2): the shared anchor rejects `delta` 0.048, `angle` 0.050, `shape` 0.042 (identical in both profiles — the anchor cell is the same cell); translation at the profile's largest effect rejects 0.040 / 0.038 / 0.028 at e = 0.1 (magnitude profile) and 0.056 / 0.056 / 0.060 at e = 1.0 (geometry profile). The two-sided entry flags the magnitude profile's translation `shape` 0.028 as not met — a conservative rate two empirical SEs below α, the known property of the two-sided form, not an excess; every one of these rates is under the upper bound.

## 6. Reproduction checks

Stated in §2: the sibling profiles reproduce each other's anchors at every column with identical p-values (3,000 / 3,000), and each slice profile reproduces the pilot's datasets at replicate indices 0–99 at every column (600 / 600 per profile, seeds and observed statistics identical, p-values not comparable across 999 and 199 permutations and excluded from the verdict). Pilot and slice replicates are never pooled.

## 7. Design-guidance statement

For a two-group, three-stage design measured as the production pipeline measures it, at the core stage effect (δ = 2.0; a saturated stage signal) and with the structural-analogue blocks:

- **Size.** Every statistic holds its nominal size at every sample size and block count (anchor ≤ 0.062 against 0.0695), and the location-only control holds the `angle` bound at every column at e = 1.0. A rejection beside a location difference still warrants care at n = 80 with two blocks, where the smaller location offsets gave `angle` 0.074–0.078.
- **`delta` (magnitude)** is the component a small cohort detects. With two blocks: e = 0.1 at 0.62 (n = 80), 0.94 (n = 200) and 1.00 (n = 500); e = 0.05 at 0.23, 0.69 and 0.96; e = 0.02 at 0.04, 0.23 and 0.63. A third block raises every rate: 0.84 / 0.99 / 1.00 at e = 0.1, 0.51 / 0.88 / 0.99 at e = 0.05, 0.11 / 0.43 / 0.83 at e = 0.02.
- **`shape`** needs n or blocks. With two blocks it stays ≤ 0.20 at n = 80, reaches 0.42–0.56 at n = 200 (e = 0.5–1.0) and 0.75–0.85 at n = 500; with three blocks 0.32–0.48 at n = 80, 0.67–0.86 at n = 200 and 0.89–0.97 at n = 500. A `shape` non-rejection at n = 80 is uninformative.
- **`angle` (orientation)**, every reading at a selected rank of 2–3: 0.10–0.15 at n = 80 and 0.16–0.24 at n = 200 with two blocks (a non-rejection is uninformative there), 0.36–0.43 at n = 500 with two blocks; 0.24–0.31, 0.41–0.49 and 0.66–0.74 with three blocks. The rate does not increase with the orientation effect between e = 0.25 and 1.0; what moves it is n and the block count. Report the cohort's selected rank beside any `angle` result.
- **Blocks.** Measuring the third block raises `delta`, `shape` and `angle` power at every n — by 0.2–0.3 at the mid-range cells — at no cost in size; block count is confounded with feature count (§9).

None of this is claimed below the core stage effect; the stage-effect arm (§8) measures how the statement degrades as δ falls toward the SEA-AD regime.

## 8. Stage-effect arm

Run on the user's decision on 2026-10-29, once the geometry slice's numbers were in (jobs 930041 and 930042, 20 + 46 shards submitted at 13:42 with the overloaded node excluded; both arrays `COMPLETED 0:0`; 7,000 and 16,000 records; 0 duplicates, missing records, signature mismatches, failures or censored surgeries; 68.1 recorded core-hours at a median 10.5–10.7 s per unit; `PROVENANCE.txt` per profile). The arm is the slice profile at n = 200, three stages and two blocks with the joint per-omic stage effect `delta_methyl = delta_expr = delta_protein` at δ ∈ {0.25, 0.5, 1.0} as explicit design points. δ = 0.25 is its configuration baseline — so no column repeats the slice, whose n = 200 / two-block column is the δ = 2.0 reference in the table — and its Type I baseline phase is a fresh 500-replicate draw at δ = 0.25.

| δ | CV stage AUROC (anchor; orientation e = 1.0) | selected rank, median (range): anchor; orientation e = 1.0 | `delta` e 0.02 / 0.05 / 0.1 | `angle` e 0.25 / 0.5 / 1.0 | `shape` e 0.25 / 0.5 / 1.0 | anchor `delta` / `angle` / `shape` | translation `angle` e = 1.0 |
|---|---|---|---|---|---|---|---|
| 0.25 | 0.997; 0.950 | 11 (6–19); 16 (5–19) | 0.050 / 0.056 / 0.084 | 0.064 / 0.088 / 0.126 | 0.048 / 0.054 / 0.046 | 0.052 / 0.058 / 0.062 | 0.054 holds |
| 0.5 | 1.000; 1.000 | 4 (3–6); 9 (6–19) | 0.068 / 0.110 / 0.310 | 0.074 / 0.110 / 0.462 | 0.036 / 0.052 / 0.066 | 0.046 / 0.050 / 0.038 | 0.048 holds |
| 1.0 | 1.000; 1.000 | 3 (2–3); 4 (3–5) | 0.058 / 0.196 / 0.508 | 0.074 / 0.070 / 0.470 | 0.070 / 0.100 / 0.168 | 0.028 / 0.048 / 0.028 | 0.048 holds |
| 2.0 (slice, §3) | 1.000; 1.000 | 2 (2–3); 3 (2–3) | 0.228 / 0.694 / 0.936 | 0.158 / 0.242 / 0.184 | 0.236 / 0.418 / 0.562 | 0.040 / 0.060 / 0.034 | 0.052 holds |

Rates ± ≤ 0.022 (R = 500); the full tables with every cell's AUROC and rank are in `tables.md`. The median rank at the `angle` cells is 12 / 13 / 16 at δ = 0.25, 4 / 5 / 9 at δ = 0.5 and 3 / 3 / 4 at δ = 1.0 for e = 0.25 / 0.5 / 1.0.

- **Size holds at every δ.** Anchor ≤ 0.062 (the δ = 0.25 `shape`) against 0.0695 and `size_holds` at all three points; translation's `angle` at e = 1.0 meets the column-scoped rule at 0.048–0.054; no translation rate at any effect exceeds 0.062. The Type I baseline phase at δ = 0.25 gives anchor 0.026 / 0.062 / 0.036 and translation 0.036 / 0.056 / 0.034 (e = 0.1, magnitude profile) and 0.028 / 0.054 / 0.060 (e = 1.0, geometry profile) — all under the upper bound; the two-sided entry flags the conservative 0.026 and 0.028 rates, not an excess.
- **The stage signal stays strong down to δ = 0.25** (held-out CV AUROC 0.997, against 0.64–0.70 in the SEA-AD cohorts and 0.55 / 0.73 at the pilot's δ = 0.05 / 0.1), so the arm covers the upper part of the transition; the orientation surgery at e = 1.0 is the only cell that erodes it (0.950 at δ = 0.25).
- **The CV rank inflates as δ falls**: a median 11 components at δ = 0.25 (up to 19), 4 at δ = 0.5, 3 at δ = 1.0 and 2 at δ = 2.0; at the orientation e = 1.0 cells 16, 9, 4 and 3.
- **`delta` and `shape` lose their power with δ.** At e = 0.1 `delta` falls from 0.94 (δ = 2.0) to 0.51, 0.31 and the null (0.08); e ≤ 0.05 is at the null for δ ≤ 0.5 and 0.20 at δ = 1.0 and e = 0.05. `shape` is at the null for δ ≤ 0.5 and reaches 0.17 at δ = 1.0 and e = 1.0 against 0.56 at δ = 2.0.
- **`angle` does not follow δ monotonically.** At e ≤ 0.5 it sits at the null for every δ ≤ 1.0 (0.06–0.11), but at e = 1.0 it reads 0.46 and 0.47 at δ = 0.5 and 1.0 — above the 0.18 at δ = 2.0 — with ranks 9 and 4 against 3. This is the rank confound of §9 in the other direction: a weaker stage signal makes the CV retain more components, which gives the largest orientation surgery room to be resolved, while a saturated signal collapses the latent space to two or three stage-driven axes. At δ = 0.25 (rank 16) the e = 1.0 rate is 0.13.
- **Against the pilot arms** (100 × 199) the arm's rates are within 0.06 (`delta`), 0.04 (`angle`) and 0.03 (`shape`) at every cell and the median ranks agree; the arm supersedes them at these three columns.

**Reproduction.** The two arm profiles' anchors reproduce each other at all three columns, 1,500 / 1,500 pairs with identical p-values (`reproduction/siblings_stage_effect.csv`). Against the pilot arms profiles the committed `--by-column` tool pairs nothing (`pilot_stage_effect_*.csv`, 0 pairs): it keys a column on its full design point, and the pilot arms' points carry seven axes (the other one-factor arms at their baseline values) against the arm's three. The seed identity holds nonetheless — the profile tests assert it at enumeration for every arm cell (replicates 0–99) — and a projected join of the records on (δ, replicate index), restricted to the pilot columns whose other axes sit at the baseline, confirms it on the data: 300 / 300 pairs per profile (100 at each δ), with generator seeds, selected ranks and observed `delta`/`angle`/`shape` identical on every pair (`reproduction/pilot_stage_effect_*_projected.csv`; p-values not comparable across 999 and 199 permutations). An option to project the column key onto named axes is the natural follow-up for the tool.

**Guidance below the core stage effect** (n = 200, three stages, two blocks): at δ = 0.25 no component has power at any effect in the enumerated ranges except a weak e = 1.0 orientation reading (0.13 at rank 16); at δ = 0.5–1.0 only `delta` at e = 0.1 (0.31–0.51) and `angle` at e = 1.0 (0.46–0.47, at ranks 9 and 4) leave the null; `shape` needs the core stage effect. A real design's held-out stage AUROC places it on this axis: 1.000 reads as δ ≥ 0.5, 0.997 as δ ≈ 0.25, and the SEA-AD level (0.64–0.70) sits below the arm, where the pilot measured no power for any component at δ ≤ 0.1.

## 9. Limitations

- **Translation at three effects.** A configuration has one effect list, so translation ran at e ∈ {0.25, 0.5, 1.0}; the predeclared rule judges e = 1.0 only, and the two smaller offsets are descriptive controls (above the bound for `angle` at n = 80 with two blocks). A per-mode effect list in the engine is the documented alternative.
- **Block count is confounded with feature count** (498 measured features with two blocks against 498 plus the protein block with three): the block effect is a measurement contrast on the same datasets, not a claim about protein data.
- **Structural analogue.** The blocks are InterSIM-calibrated structural analogues of methylation, expression and protein data; the rates transfer to a real design only through the stage signal and the feature counts, not through any biology.
- **Saturated stage signal.** At δ = 2.0 the held-out CV stage AUROC is 1.000 everywhere, far above the SEA-AD cohorts (0.64–0.70); the statement in §7 holds for a strongly staged design, and §8 is the only part of this report that speaks to weaker staging.
- **The rank confound of orientation power.** `angle` is measured in the CV-selected latent space (rank 2–3 here) and its power follows n and the block count rather than the orientation effect; a flat or non-monotone response in e is a property of the measurement, not evidence that larger orientation differences are harder to detect. The latent-rank ladder (2026-09-08) and the pilot report §4 document the mechanism; every `angle` reading above carries its rank.
- **Scope.** Three stages, ρ = 0, equal groups, `p_dmp` 0.1, n ≤ 500. The pilot's four-stage, n = 1000, group-ratio, stage-proportion, `p_dmp` and ρ columns stand at pilot grade (100 × 199) and are not rerun here. Group-specific dispersion is not in the generator (change B of the pilot's follow-ups).
