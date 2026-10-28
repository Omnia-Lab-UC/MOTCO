# Design-guidance pilot — findings report (2026-10-28)

**Change:** `design-guidance-grid` · **Configs:** `examples/trajectory_power_study/design_guidance_{core,arms}_{magnitude,geometry}_pilot.json` · **Code revision:** `7c218c6` (branch `feat/design-guidance-grid`) · **Records:** `results/design-guidance-pilot-2026-10-28/<profile>/` (`report/`, `PROVENANCE.txt`) and `results/design-guidance-pilot-2026-10-28/figures/` · **Grade:** pilot, 100 replicates × 199 permutations; every rate below carries its Monte Carlo standard error at R = 100 (SE ≤ 0.05) and the anchor bound is α + 2·√(α(1 − α)/R) = 0.0936. **This report proposes the slice to rerun at paper grade; it makes no paper-grade claim.**

## 1. Configuration and provenance

_Pending the run: cite each profile's `PROVENANCE.txt` (date, config sha256, code revision, host and sbatch line with `STUDY_N_JOBS` unset, BLAS threads, units, versions, job ids, outcome, wall and unit timings, shard layout, error policy, merge and report commands)._

Shared design (one baseline, one seed root `base_seed` 1028, one matched-seed family `design-guidance`): n = 200, three stages, `group_ratio` 0.5, `p_dmp` 0.1, ρ = 0, every per-omic δ 2.0, `magnitude_kind: joint`; pooled PLS on methylation + expression with the stage-supervised double-CV rank (`cv1_splits` 3, `cv2_splits` 4, 5 repeats, ≤ 20 components, `random_state` 1203), `orientation_estimator: net_displacement`, `min_components: 2`, 199 RRPP permutations, `n_jobs` 1. The production measurement is named in every profile, not defaulted. Attribution and the Phase 4 gate are off; the only predeclared rule is the advisory `anchor_decision` (α 0.05, k = 2, `delta`/`angle`/`shape`).

- **Core** (design grid `n_samples` {80, 200, 500, 1000} × `n_stages` {3, 4} × `layers` {methylation + expression, all three}): magnitude at e ∈ {0.01, 0.02, 0.05, 0.1, 0.25} (98 cells, 9,800 units); orientation, shape and translation at e ∈ {0.25, 0.5, 1.0} (162 cells, 16,200 units).
- **Arms** (explicit points, one setting moved off the baseline): group ratio {0.65, 0.8}; stage proportions 0.5/0.3/0.2; `p_dmp` {0.05, 0.15}; joint δ {0.05, 0.1, 0.25, 0.5, 1.0}; ρ {0.5, 0.8}. Magnitude 80 cells / 8,000 units; geometry 132 cells / 13,200 units.

The magnitude and geometry siblings share their seed family, so at every design column their zero-effect anchors evaluate the same datasets with the same RRPP draws (the smoke verified identical p-values at every column). This report reads the **geometry profile's anchor** and uses the magnitude profile's only as the reproduction check in §4; the two are never pooled.

## 2. Unit and failure accounting

_Pending the run: units expected vs present once, signature mismatches, failed records per column (especially n = 80 with four stages), censored surgeries (`realized_surgery.csv`)._

## 3. Per-column operating characteristics

_Pending the run: the tables from `report/design_point_operating.csv` (target-statistic power per column and effect, rate ± MC SE), the four figures from `scripts/design_guidance_figure.py`, and the reading per component (magnitude/`delta`, orientation/`angle`, shape/`shape`) across sample size, stage count and block count, then across the arms._

## 4. Stage signal and selected rank

_Pending the run: each column's held-out CV stage AUROC (`report/design_point_stage_signal.csv`) and selected rank (`median_selected_lv`), read against the SEA-AD cohorts' recorded 0.64–0.70. The smoke already showed AUROC 1.000 at δ = 2.0 at every core column and 0.988 at δ = 0.25, which is why the δ arm was extended to 0.05 and 0.1 at the approval gate._

## 5. Size at every column

_Pending the run: every column's zero-effect anchor rate per statistic with the anchor-decision status (`report/anchor_decision.json`), translation's `angle` rate at e = 1.0 read against the 0.0936 bound (it sat at 0.058–0.068 against the paper-grade bound 0.0695 in the net-displacement studies), and the magnitude profile's reproduction of the geometry profile's anchor._

## 6. The originally planned grid and what was run instead

The project's originally planned simulation grid (Table 1) crossed every row below; crossing them is infeasible (3,456 null combinations before any power cell), its mean shifts (0.01–0.2) sit far below the InterSIM stage effect the generator operates on (δ = 2.0), its "power" rows do not map to trajectory components, and two of them are null scenarios. Each factor's disposition:

| Table 1 factor | Disposition | Where, and why |
|---|---|---|
| Sample size | **Run** | Core axis `n_samples` ∈ {80, 200, 500, 1000}, crossed with stages and blocks. |
| Group ratio | **Run** | Arm: `group_ratio` 0.65 and 0.8 off 0.5. |
| Number of stages | **Run** | Core axis `n_stages` ∈ {3, 4}. Two-stage designs are dropped: with two stage means the configuration is a segment, `shape` is undefined by construction (nothing to bend) and `angle` reduces to the segment direction, so the three-component guidance cannot be read there; the case study has three stages. |
| Stage proportions | **Run** | Arm: 0.5/0.3/0.2 off equal thirds. |
| Differential-feature proportion | **Run** | Arm: `p_dmp` 0.05 and 0.15 off 0.1. 0.2 is excluded because at three stages it fails the surgery-headroom check (orientation saturates at e = 0.89, translation at 0.52). |
| Shared stage shift | **Run, rescaled** | Arm: the joint per-omic δ ∈ {0.05, 0.1, 0.25, 0.5, 1.0} off InterSIM's 2.0. Table 1's 0.01–0.2 are on a scale the generator never reaches; 0.05–0.1 is where the rank-floor check found the SEA-AD stage-signal regime. |
| Noise level | **Dropped** | Redundant with the stage shift: features are standardized before the PLS fit, so σ and δ enter the measurement only through δ/σ, which the δ arm already sweeps. |
| Group-specific mean shift | **Run** | The `translation` mode (a constant location offset of group B) at e ∈ {0.25, 0.5, 1.0} in both geometry profiles, read as the negative control for all three statistics; the `magnitude` mode covers a group-specific *size* change. |
| Covariance scaling | **Deferred** | To the group-specific-dispersion change (change B): an opt-in generator field scaling group B's covariance and a null-only Type I study crossing spread {1, 1.5, 2, 3} × group ratio {0.2, 0.5, 0.8}. It is a size question, not a power row. |
| Correlation scaling | **Deferred** | Same change and the same reason: it changes group B's dispersion, not its trajectory. |
| Feature filtering | **Dropped** | Supervised filtering (by group or stage) biases the permutation null unless it is repeated inside every permutation, so MOTCO documents it as forbidden; label-blind filtering is upstream preprocessing and changes what is measured, not the validity of the test. |
| Integration method | **Dropped** | PLS is the production measurement (`docs/roadmap.md`); SNF with the Euclidean statistics is deferred pending graph-native statistics; MOFA is out of scope; `concat` is a diagnostic baseline, not a constructed latent space. |

## 7. Proposed paper-grade slice

_Pending the run._

## 8. Limitations

- **Pilot precision.** 100 replicates: SE ≤ 0.05 on every rate, and the anchor bound (0.0936) is loose; a column whose anchor sits between 0.07 and 0.09 is neither cleared nor flagged at this grade. Translation's `angle` rate is read against the bound and flagged, never claimed.
- **Block count is confounded with feature count.** The two-block columns measure 498 features (methylation + expression); the three-block columns add the protein block. The contrast is "more blocks and more features", not block count alone; block weighting is roadmap item 5.
- **Structural analogue.** The InterSIM blocks are a structural analogue of a real multi-omic design (498 measured features, cross-omic coupling through the reference maps), not an emulation of SEA-AD's RNA/ATAC.
- **Stage signal.** The core columns run at δ = 2.0, where the held-out stage AUROC is saturated; only the δ arm reaches the regime of a real analysis, and only at n = 200 with three stages and two blocks.
- **One-at-a-time arms** measure main effects off one baseline; interactions with sample size are not measured.
