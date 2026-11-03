# Presentation figures: speaker notes

Rendered 2026-11-02T20:23:34-03:00 at 849534f82550efbb92caac137540ffb7173cf4f3 (clean worktree), MOTCO 0.14.0, matplotlib 3.10.9, seed 0. Regenerate with `python scripts/presentation_figures.py --out-dir <dir>`; `manifest.json` holds the same content with every source's SHA-256.

## 1. What each trajectory mode does (`mode_geometry`)

**Message.** Each simulated mode changes group B's trajectory in one geometric way, read by one statistic: magnitude scales it (size δ), orientation rotates it (θ), shape bends it (shape), and translation moves it without changing any of the three.

**Measurement.** None: a constructed schematic, no data and no source file.

**Caveats.**

- Constructed coordinates, not data or a fitted space.
- Magnitude scales by 1.6 and orientation rotates by 35° about stage 1; shape moves the interior stage only, so it needs at least three stages.

**Files.** `mode_geometry.png`, `mode_geometry.pdf`

## 2. The four modes as the generator realizes them (`mode_simulated`)

**Message.** The numpy generator realizes each mode on a shared baseline: group B's stage-mean path is scaled, rotated, bent or shifted against group A's, and coincides with it at zero effect.

**Measurement and sources.**

- generator input (no measurement) (`src/motco/simulations/data/intersim_reference.npz`, numpy generator reference data)

**Caveats.**

- Display projection: per panel, the first two principal axes of the six group × stage means of the concatenated standardized blocks; translation's offset is orthogonal to the stage geometry, so its panel's second axis is the group offset. It is not the measurement space (the double-CV PLS latent space).
- One dataset per mode from generate_showcase_datasets (seed given in the manifest, n = 1,000, 3 stages, e = 1.0, surgery_censoring='clamp' as in the showcase); the modes share group A's baseline draw but are not otherwise identical.
- Axes are flipped so group A runs from stage 1 toward positive values (translation: B above A); the sign carries nothing.
- Rendered with seed 0 and n = 1,000 per dataset.

**Files.** `mode_simulated.png`, `mode_simulated.pdf`, `mode_simulated_anchor.png`, `mode_simulated_anchor.pdf`, `mode_simulated_magnitude.png`, `mode_simulated_magnitude.pdf`, `mode_simulated_orientation.png`, `mode_simulated_orientation.pdf`, `mode_simulated_shape.png`, `mode_simulated_shape.pdf`, `mode_simulated_translation.png`, `mode_simulated_translation.pdf`

## 3. Size holds from n = 80 to n = 1,501 (`size_atlas`)

**Message.** Every zero-effect anchor and translation control at every committed design, from n = 80 to n = 1,501 and with up to six covariate columns, rejects at or below its own α + k·SE bound.

**Measurement and sources.**

- PLS, CV rank ≥ 2 · net displacement · joint magnitude · uncorrected RRPP · 500 × 999 (`results/design-guidance-slice-2026-10-29/core_geometry/report/anchor_decision.csv`, design-guidance-slice-2026-10-29 · core_geometry)
- PLS, CV rank ≥ 2 · net displacement · joint magnitude · uncorrected RRPP · 500 × 999 (`results/design-guidance-slice-2026-10-29/stage_effect_geometry/report/anchor_decision.csv`, design-guidance-slice-2026-10-29 · stage_effect_geometry)
- PLS, CV rank ≥ 2 by floor check (0/152 ranks change) · net displacement · joint magnitude · uncorrected RRPP · 500 × 999 (`results/net-displacement-type-i-2026-10-13/phase5/report/anchor_decision.csv`, net-displacement-type-i-2026-10-13 · phase5)
- PLS, CV rank ≥ 2 · net displacement · joint magnitude · uncorrected RRPP · 500 × 999 (`results/rank-floor-2026-10-23/delta0.1/rank_floor_rates.csv`, rank-floor-2026-10-23 · delta0.1 (floor_2 arm))
- PLS, CV rank ≥ 2 · net displacement · joint magnitude · uncorrected RRPP · 500 × 999 (`results/rank-floor-2026-10-23/delta0.05/rank_floor_rates.csv`, rank-floor-2026-10-23 · delta0.05 (floor_2 arm))
- PLS, CV rank ≥ 2 · net displacement · joint magnitude · 6 covariate columns · df-corrected RRPP · 500 × 999 (`results/corrected-permutation-2026-10-30/size_readings.csv`, corrected-permutation-2026-10-30 (df arm))
- PLS, CV rank ≥ 2 · net displacement · joint magnitude · 6 covariate columns · uncorrected RRPP · 500 × 999 (`results/adni-design-2026-11-02/run/size_readings.csv`, adni-design-2026-11-02 · run)
- latent-space fits at floor 1 and floor 2 on rebuilt datasets; no RRPP (`results/rank-floor-equivalence-2026-11-02/summary.json`, rank-floor-equivalence-2026-11-02)

**Caveats.**

- Bounds are those the sources record: α + 2·SE at R = 500 (0.0695) except the corrected-permutation readings, judged per reading at α + 3·SE (0.0792).
- The harness tables record no verdict for translation; its bound is the anchor's form at its own R and its verdict is derived against it.
- Production arm per source: floor 2 (rank-floor check), df (corrected permutation), adjusted (ADNI), the net-displacement CV column (n = 1,200).
- Excluded: case-study-type-i (superseded by the df configuration) and covariate-type-i (pre-floor), whose uncorrected excess is what covariates_null shows; the latent-adjustment arms (a sensitivity analysis); the PC1 columns; the Phase 6 Type I profile (pre-floor, no equivalence check; the same table at the current measurement appears through the rank-floor and corrected-permutation rows).

**Files.** `size_atlas.png`, `size_atlas.pdf`

## 4. At n = 1,200 every statistic detects its mode (`specificity_n1200`)

**Message.** At n = 1,200 each mode is detected by its target statistic; magnitude moves size alone, while rotating or bending a trajectory also changes its size (and rotation its shape).

**Measurement and sources.**

- PLS, CV rank ≥ 2 by floor check (0/152 ranks change) · net displacement · joint magnitude · uncorrected RRPP · 500 × 999 (`results/net-displacement-power-2026-10-14/phase5/report/design_point_operating.csv`, net-displacement-power-2026-10-14 · phase5)
- latent-space fits at floor 1 and floor 2 on rebuilt datasets; no RRPP (`results/rank-floor-equivalence-2026-11-02/summary.json`, rank-floor-equivalence-2026-11-02)

**Caveats.**

- Phase 5 net-displacement power run (2026-10-14), net-displacement column; rows at e = 1.0, the anchor at e = 0.
- Magnitude e = 1.0 here is not comparable with e = 0.1 at n = 80; the message is the off-diagonal pattern.

**Files.** `specificity_n1200.png`, `specificity_n1200.pdf`

## 5. At n = 80 the statistics mix (`specificity_n80`)

**Message.** At n = 80 (the SEA-AD scale) only magnitude stays specific: orientation moves all three statistics about equally and shape moves size more than shape, so one significant statistic does not name the mode.

**Measurement and sources.**

- PLS, CV rank ≥ 2 · net displacement · joint magnitude · uncorrected RRPP · 500 × 999 (`results/design-guidance-slice-2026-10-29/core_magnitude/report/design_point_operating.csv`, design-guidance-slice-2026-10-29 · core_magnitude)
- PLS, CV rank ≥ 2 · net displacement · joint magnitude · uncorrected RRPP · 500 × 999 (`results/design-guidance-slice-2026-10-29/core_geometry/report/design_point_operating.csv`, design-guidance-slice-2026-10-29 · core_geometry)

**Caveats.**

- Design-guidance slice core profiles at n = 80, 3 stages, methylation + expression; magnitude at e = 0.1, orientation, shape and translation at e = 1.0, the anchor at e = 0.
- Same colour scale (0 to 1) as specificity_n1200.

**Files.** `specificity_n80.png`, `specificity_n80.pdf`

## 6. Power against sample size (`power_vs_n`)

**Message.** Power rises with n for every statistic and with the third block: size (δ) is well powered by n = 200 at e = 0.1, shape needs n ≈ 500 with two blocks, and orientation stays below 0.75 even at n = 500.

**Measurement and sources.**

- PLS, CV rank ≥ 2 · net displacement · joint magnitude · uncorrected RRPP · 500 × 999 (`results/design-guidance-slice-2026-10-29/figures/design_guidance_power.csv`, design-guidance-slice-2026-10-29 · figures)

**Caveats.**

- Design-guidance slice (2026-10-29), 3 stages, 500 × 999 per cell; ±1 Monte Carlo SE.
- Orientation responds little to e at this rank (2–3); see the detail figure.

**Files.** `power_vs_n.png`, `power_vs_n.pdf`, `power_vs_n_delta.png`, `power_vs_n_delta.pdf`, `power_vs_n_angle.png`, `power_vs_n_angle.pdf`, `power_vs_n_shape.png`, `power_vs_n_shape.pdf`

## 7. A weak stage signal inflates the rank (`power_vs_stage_signal`)

**Message.** As the joint stage effect weakens from δ = 2 to 0.25, the CV-selected rank inflates (median 2 to 11) and size and shape lose their power while the stage stays almost perfectly classified; orientation peaks at intermediate δ.

**Measurement and sources.**

- PLS, CV rank ≥ 2 · net displacement · joint magnitude · uncorrected RRPP · 500 × 999 (`results/design-guidance-slice-2026-10-29/stage_effect_magnitude/report/design_point_operating.csv`, design-guidance-slice-2026-10-29 · stage_effect_magnitude)
- PLS, CV rank ≥ 2 · net displacement · joint magnitude · uncorrected RRPP · 500 × 999 (`results/design-guidance-slice-2026-10-29/stage_effect_magnitude/report/design_point_stage_signal.csv`, design-guidance-slice-2026-10-29 · stage_effect_magnitude)
- PLS, CV rank ≥ 2 · net displacement · joint magnitude · uncorrected RRPP · 500 × 999 (`results/design-guidance-slice-2026-10-29/stage_effect_geometry/report/design_point_operating.csv`, design-guidance-slice-2026-10-29 · stage_effect_geometry)
- PLS, CV rank ≥ 2 · net displacement · joint magnitude · uncorrected RRPP · 500 × 999 (`results/design-guidance-slice-2026-10-29/stage_effect_geometry/report/design_point_stage_signal.csv`, design-guidance-slice-2026-10-29 · stage_effect_geometry)
- PLS, CV rank ≥ 2 · net displacement · joint magnitude · uncorrected RRPP · 500 × 999 (`results/design-guidance-slice-2026-10-29/core_magnitude/report/design_point_operating.csv`, design-guidance-slice-2026-10-29 · core_magnitude)
- PLS, CV rank ≥ 2 · net displacement · joint magnitude · uncorrected RRPP · 500 × 999 (`results/design-guidance-slice-2026-10-29/core_magnitude/report/design_point_stage_signal.csv`, design-guidance-slice-2026-10-29 · core_magnitude)
- PLS, CV rank ≥ 2 · net displacement · joint magnitude · uncorrected RRPP · 500 × 999 (`results/design-guidance-slice-2026-10-29/core_geometry/report/design_point_operating.csv`, design-guidance-slice-2026-10-29 · core_geometry)
- PLS, CV rank ≥ 2 · net displacement · joint magnitude · uncorrected RRPP · 500 × 999 (`results/design-guidance-slice-2026-10-29/core_geometry/report/design_point_stage_signal.csv`, design-guidance-slice-2026-10-29 · core_geometry)

**Caveats.**

- n = 200, 3 stages, two blocks; δ = 0.25–1.0 from the slice's stage-effect arm, δ = 2 from the slice core at the same design.
- Stage AUROC is the median held-out CV AUROC of the rank-selection fits; it stays near 1 in this range (the pilot's δ = 0.05 and 0.1 reach 0.55–0.73, bracketing SEA-AD).

**Files.** `power_vs_stage_signal.png`, `power_vs_stage_signal.pdf`

## 8. Net displacement measures orientation better (`estimator_comparison`)

**Message.** On identical datasets, net displacement doubles orientation power at n = 80 (0.300 → 0.648 at e = 1) and lifts it at n = 1,200 (0.840 → 0.924), while removing the shape mode's leak into θ.

**Measurement and sources.**

- n = 80 panel: PLS, CV rank ≥ 1 (pre-floor) · PC1 and net displacement on the same datasets · joint magnitude · uncorrected RRPP · 500 × 999 (`results/net-displacement-power-2026-10-14/phase6/report/design_point_operating.csv`, net-displacement-power-2026-10-14 · phase6)
- n = 1,200 panel: PLS, CV rank ≥ 2 by floor check (0/152 ranks change) · PC1 and net displacement on the same datasets · joint magnitude · uncorrected RRPP · 500 × 999 (`results/net-displacement-power-2026-10-14/phase5/report/design_point_operating.csv`, net-displacement-power-2026-10-14 · phase5)
- latent-space fits at floor 1 and floor 2 on rebuilt datasets; no RRPP (`results/rank-floor-equivalence-2026-11-02/summary.json`, rank-floor-equivalence-2026-11-02)

**Caveats.**

- Both designs predate the rank floor; the n = 1,200 panel carries the floor check's tag, the n = 80 panel stays pre-floor (the comparison is between estimators, not a current-measurement claim).
- Points at e = 0 are the shared zero-effect anchor.

**Files.** `estimator_comparison.png`, `estimator_comparison.pdf`

## 9. Covariates and the permutation null (`covariates_null`)

**Message.** Each covariate column shrinks the residual df and inflates the uncorrected δ null (pooled 0.0545 at one column, 0.066 at six, above the 0.0597 bound); the df-corrected permutation holds it (0.0545 at six).

**Measurement and sources.**

- PLS, CV rank ≥ 2 · net displacement · joint magnitude · 1 / 3 / 6 covariate columns · uncorrected vs df-corrected RRPP · 500 × 999 (`results/case-study-diagnosis-2026-10-30/delta0.1/diagnosis_rates.csv`, case-study-diagnosis-2026-10-30 · delta0.1)
- PLS, CV rank ≥ 2 · net displacement · joint magnitude · 1 / 3 / 6 covariate columns · uncorrected vs df-corrected RRPP · 500 × 999 (`results/case-study-diagnosis-2026-10-30/delta0.05/diagnosis_rates.csv`, case-study-diagnosis-2026-10-30 · delta0.05)

**Caveats.**

- Case-study diagnosis (2026-10-30), covariate-blind basis (z_noise arms), δ = 0.1 and 0.05 × anchor and translation; pooled over 2,000 replicates per point.
- Faint points are the four individual readings (R = 500 each, bound 0.0695).

**Files.** `covariates_null.png`, `covariates_null.pdf`
