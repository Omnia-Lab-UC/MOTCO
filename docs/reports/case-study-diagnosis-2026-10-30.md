# Case-study configuration diagnosis — 2026-10-30

> This report diagnoses the size excess that the [case-study configuration Type I check](case-study-type-i-2026-10-29.md)
> found at the SEA-AD case study's configuration (openspec change `case-study-configuration-diagnosis`, spec
> `case-study-configuration-diagnosis`). It covers two findings: the covariate-adjusted model's `delta`, and the
> `stage` latent adjustment's anchor `angle` and `shape`.
>
> Every number below is read from `results/case-study-diagnosis-2026-10-30/`:
> - `case_study_diagnosis_decision.json`, `diagnosis_contrasts.csv`, `decide.json`;
> - per run in `delta0.1/` and `delta0.05/`: `diagnosis_rates.csv`, `diagnosis_rank_table.csv`,
>   `diagnosis_reproduction.json`, `merge.json`, `settings.json`, `PROVENANCE.txt`;
> - the smoke records `smoke-delta0.1/`, `smoke-delta0.05/` and `smoke-decision/`.
>
> `diagnosis_pvalues.csv`, `merged.jsonl` and `records/` stay on `ing`; each `PROVENANCE.txt` names them with their
> sha256. Rates are rejection rates at p ≤ 0.05 over R = 500 replicates at 999 permutations. The bound is the parent
> check's one-sided `α + 2·√(α(1 − α)/R)` = 0.0695, i.e. at most 34 rejections in 500.

**Run.**
- Two 50-shard SLURM arrays on `ing` (partition `512x1024`, n11 and the heavily loaded n3 excluded), from the clean
  worktree of `feat/case-study-configuration-diagnosis` at `5c30c51`. The git tree is `4fe2050d…`, identical to the
  workstation, and the diagnosis tests passed there.
- Job 930358 ran the δ = 0.1 parent and job 930408 the δ = 0.05 parent: 50 / 50 tasks `COMPLETED 0:0` in each,
  10.7 + 14.4 = 25.1 recorded core-hours (34 planned).
- 1,000 units per parent (the zero-effect anchor and the translation control at e = 1.0, × 500 replicates).
  Each unit ran four analyses and twelve RRPP-only arms, with 0 arm failures.
- Every unit passed three hard checks:
  - the reconstruction of the fitted covariate-blind projection matched the analysis's latent scores (maximum
    error 2.7e-13 against the 1e-9 tolerance);
  - `X − X0` and `Y − Y0` were rank one in `z` (maximum residual 6.5e-16);
  - the script's `z_noise5` equalled `analyze_trajectories` exactly.

**Reproduction.** The diagnosis regenerated the parent's data, covariates, CV folds and RRPP seeds. `z_noise5`
reproduces the parent's `adjusted` arm, and `latent_stage_cv` its `latent_stage` arm, with **0 mismatches over 3,000 +
3,000 p-values per parent** (`diagnosis_reproduction.json`, status `holds` in both runs). So every derived arm below sits
on the parent's exact units, and the decision's gate is open.

**Verdicts.** Both are read verbatim from `case_study_diagnosis_decision.json` (`"decision": "diagnosed"`,
`"indicative": false`):
- **Covariate-adjusted `delta`: `residual_df_shrinkage`.** The column contrast C2 (`z_noise5` vs `z_noise0`) is
  significant (one-sided p = 0.015). The residual-df corrected case-study model holds the bound for every statistic at
  both controls in both runs, and its own column contrast C4 is not significant (p = 0.156).
- **`stage` latent adjustment: `residualized_basis`.** At a forced rank of 2, with no cross-validation, it still
  exceeds on anchor `angle` in both runs (0.084 / 0.082), although its `shape` now holds.

**In one paragraph.**
- **The design.** The diagnosis reran the parent check's exact units:
  - on the covariate-blind basis, it crossed the confounded covariate `z` (estimated in the model, or its injection
    removed exactly from the latent scores) with 0, 2 or 5 pure-noise covariates;
  - each of those six arms was rerun with the Freedman–Lane residuals rescaled by κ = √(n / (n − p_red));
  - the `stage` latent adjustment ran at the CV floor of 2, at floor 1 and at a forced rank of 2.
- **`delta` rises with the column count and the correction removes it.** Pooled over the four control readings
  (2,000 units), `delta` rises 0.0545 → 0.0595 → 0.066 with `z` in the model and 0.0525 → 0.0545 → 0.060 under the
  oracle. The corrected arms read 0.0475–0.0545, with no arm above the bound in any reading.
- **Estimating `z` costs no `delta` size.** C3 (`z_noise0` vs `oracle_noise0`) has p = 0.39.
- **What this means.** The case study's `delta` excess is the Freedman–Lane null's residual-df shrinkage at ten
  reduced-model columns and n = 80. It is not the covariate's confounding.
- **The `stage` arm.** Its anchor `angle` excess sits at rank 2 itself. Its `shape` excess sits in the CV rank's
  tail, and a fixed rank of 2 cures `shape` but not `angle`.
- **Consequences.** The case study's primary waits for a MOTCO change offering the corrected permutation, validated
  by its own Type I and power check. The `stage` sensitivity arm stays inadmissible at this configuration, at a
  fixed rank too.

## 1. Construction

`scripts/case_study_diagnosis.py` (design D1) binds each run to one committed run directory of the parent check
(`--parent`).

**Parent binding.**
- It refuses a parent whose `settings.json` fails its recorded hash.
- It also refuses one whose settings `case_study_type_i_check.build_settings` does not rebuild from the recorded
  arguments to the same hash.
- The diagnosis's own `settings.json` carries the parent's settings, path and hash: `ca57cc38…` against parent
  `899ce29a…` at δ = 0.1, and `6c5d447a…` against `97eb6d90…` at δ = 0.05.

The parent's construction is unchanged:
- the Phase 6 SEA-AD design, `[[11, 10, 28], [9, 10, 12]]` with n = 80, methylation + expression;
- the per-omic stage effect weakened to δ = 0.1 or 0.05;
- `z = 0.5·g + 1·s + 2·g·s + ε`, injected at `s_z = 10` along a fixed direction;
- five i.i.d. N(0, 1) noise covariates `x1`–`x5`;
- the rank floor of 2, the harness CV knobs and the net-displacement estimator.

**Arms** (design D2–D6). Every arm of a unit shares its data and its RRPP seed.

| arm | latent scores | covariates in both trajectory models | p_red | κ |
|---|---|---|---|---|
| `oracle_noise0` | `Y0` | none | 4 | — |
| `oracle_noise2` | `Y0` | `x1`, `x2` | 6 | — |
| `oracle_noise5` | `Y0` | `x1`–`x5` | 9 | — |
| `z_noise0` | `Y` | `z` | 5 | — |
| `z_noise2` | `Y` | `z`, `x1`, `x2` | 7 | — |
| `z_noise5` | `Y` | `z`, `x1`–`x5` (the case study's model) | 10 | — |
| `<arm>_corrected` | as `<arm>` | as `<arm>` | as `<arm>` | √(80 / (80 − p_red)): 1.026–1.069 |
| `latent_stage_cv` | own basis, `stage` adjustment, CV at floor 2 (the parent's arm) | `z`, `x1`–`x5` | 10 | — |
| `latent_stage_floor1` | own basis, `stage` adjustment, CV at floor 1 (descriptive) | `z`, `x1`–`x5` | 10 | — |
| `latent_stage_rank2` | own basis, `stage` adjustment, forced rank 2, no CV | `z`, `x1`–`x5` | 10 | — |

- **`Y` and the factorial.** `Y` is the latent score matrix of one `analyze_trajectories` call with the six covariates
  on the covariate-blind basis. The six factorial arms are inference-only on that basis: the analysis's design,
  `_covariate_adjusted_design`, `estimate_difference` and serial `RRPP`.
- **The oracle (design D3).**
  - `X0` is D10's build at `s_z = 0` with the same seeds.
  - The fitted covariate-blind map is reconstructed from the block scalers, the PLS centring and `ddof = 1` scaling,
    and the PLS rotations, and checked against the analysis's latent scores.
  - `Y0` is `X0` through that map. It equals `Y` minus `z` times a fixed latent vector, which the tests check against
    the analytic injection.
- **The correction (design D5).** Each corrected arm runs `RRPP` on `ŷ_red + κ·e_red` with the arm's models and seed:
  the same permutation indices, the uncorrected observed statistics, and p-values by the analysis's rule.

**Contrasts** (design D7) are paired on the same units:
- For a contrast of arm B against arm A, replicate r contributes `d_r`: the readings (2 runs × 2 controls) that B
  rejects and A does not, minus the reverse.
- The one-sided p-value is P(Σ ε_r d_r ≥ Σ d_r) over independent equiprobable signs, computed exactly by convolution.
  Each replicate is one cluster because its four readings share the covariates, noise columns and RRPP seed.
- The predeclared `delta` contrasts are:
  - C1: `oracle_noise5` vs `oracle_noise0`
  - C2: `z_noise5` vs `z_noise0`
  - C3: `z_noise0` vs `oracle_noise0`
  - C4: `z_noise5_corrected` vs `z_noise0_corrected`

**Smoke** (design D10; `smoke-*`, R = 20 × 99 permutations; the `design.md` "Smoke and shard plan" record):
- 40 / 40 units per parent with 0 failures, and indicative verdicts `unresolved` / `residualized_basis`;
- a separate R = 4 × 999 run reproduced 24 + 24 p-values per parent with 0 mismatches;
- one unit took 20.5 s on a workstation core (1.11× a parent unit). On `ing` the mean unit took 31 s (δ = 0.1,
  n2 and n4) and 42 s (δ = 0.05, n2).

## 2. Rates against the bound

Rejections out of 500 per reading, then the rate pooled over the four readings, from each run's
`diagnosis_rates.csv`. The readings are the anchor and translation controls at δ = 0.1 and 0.05. **Bold** marks a
reading above 34 / 500, i.e. verdict `exceeds`.

### `delta`

| arm | anchor δ = 0.1 | anchor δ = 0.05 | translation δ = 0.1 | translation δ = 0.05 | pooled (2,000) |
|---|---|---|---|---|---|
| `oracle_noise0` | 31 | 23 | 27 | 24 | 0.0525 |
| `oracle_noise2` | 34 | 28 | 23 | 24 | 0.0545 |
| `oracle_noise5` | **37** | 29 | 28 | 26 | 0.0600 |
| `z_noise0` | 34 | 30 | 23 | 22 | 0.0545 |
| `z_noise2` | **36** | 29 | 28 | 26 | 0.0595 |
| `z_noise5` | 34 | 32 | **36** | 30 | 0.0660 |
| `oracle_noise0_corrected` | 28 | 19 | 26 | 23 | 0.0480 |
| `oracle_noise2_corrected` | 32 | 24 | 21 | 18 | 0.0475 |
| `oracle_noise5_corrected` | 32 | 24 | 24 | 21 | 0.0505 |
| `z_noise0_corrected` | 30 | 26 | 22 | 20 | 0.0490 |
| `z_noise2_corrected` | 31 | 24 | 22 | 19 | 0.0480 |
| `z_noise5_corrected` | 30 | 24 | 28 | 27 | 0.0545 |
| `latent_stage_cv` | 30 | 32 | 24 | 26 | 0.0560 |
| `latent_stage_floor1` | 31 | 33 | 24 | 28 | 0.0580 |
| `latent_stage_rank2` | 30 | 31 | 25 | 25 | 0.0555 |

`z_noise5` reproduces the parent's `adjusted` `delta` reading for reading: 0.068 / 0.064 at the anchor and 0.072 /
0.060 at translation, including the one excess (translation δ = 0.1, 36 / 500).

### `angle` and `shape`

| arm | `angle` (four readings) | `angle` pooled | `shape` (four readings) | `shape` pooled |
|---|---|---|---|---|
| `oracle_noise0` | **40 39 40 40** | 0.0795 | 30 26 28 30 | 0.0570 |
| `oracle_noise2` | **40 39 42 43** | 0.0820 | 26 23 30 30 | 0.0545 |
| `oracle_noise5` | **41 38 43 44** | 0.0830 | 34 29 29 28 | 0.0600 |
| `z_noise0` | 26 27 24 23 | 0.0500 | 26 22 **35 35** | 0.0590 |
| `z_noise2` | 24 23 19 18 | 0.0420 | 24 25 32 29 | 0.0550 |
| `z_noise5` | 24 24 21 22 | 0.0455 | 28 27 34 34 | 0.0615 |
| `oracle_noise0_corrected` | **38 37 36 37** | 0.0740 | 27 24 27 25 | 0.0515 |
| `oracle_noise2_corrected` | **38 38 37 35** | 0.0740 | 20 18 29 26 | 0.0465 |
| `oracle_noise5_corrected` | **36 36 37 36** | 0.0725 | 27 25 27 22 | 0.0505 |
| `z_noise0_corrected` | 22 23 22 21 | 0.0440 | 23 21 32 31 | 0.0535 |
| `z_noise2_corrected` | 20 20 15 14 | 0.0345 | 18 18 28 27 | 0.0455 |
| `z_noise5_corrected` | 20 21 21 22 | 0.0420 | 19 17 28 29 | 0.0465 |
| `latent_stage_cv` | **45 42** 32 32 | 0.0755 | **46 39** 28 28 | 0.0705 |
| `latent_stage_floor1` | **41 40** 32 29 | 0.0710 | **41** 34 25 26 | 0.0630 |
| `latent_stage_rank2` | **42 41** 31 30 | 0.0720 | 34 32 27 25 | 0.0590 |

The readings are in the column order of the `delta` table. `latent_stage_cv` reproduces the parent's `latent_stage`
arm: anchor `angle` 0.090 / 0.084 and `shape` 0.092 / 0.078.

## 3. Paired contrasts

Clustered over both runs and both controls: 500 replicate clusters, 2,000 paired readings
(`diagnosis_contrasts.csv`, scope `clustered`). "B only" and "A only" count discordant readings, and p is the exact
one-sided sign-flip p-value that B rejects more often than A. The `delta` rows are judged and the rest descriptive.

| contrast | statistic | B | A | B rejections | A rejections | B only | A only | net | p |
|---|---|---|---|---|---|---|---|---|---|
| C1 | `delta` | `oracle_noise5` | `oracle_noise0` | 120 | 105 | 33 | 18 | 15 | 0.065 |
| C2 | `delta` | `z_noise5` | `z_noise0` | 132 | 109 | 41 | 18 | 23 | **0.015** |
| C3 | `delta` | `z_noise0` | `oracle_noise0` | 109 | 105 | 31 | 27 | 4 | 0.386 |
| C4 | `delta` | `z_noise5_corrected` | `z_noise0_corrected` | 109 | 98 | 32 | 21 | 11 | 0.156 |
| C1 | `angle` | `oracle_noise5` | `oracle_noise0` | 166 | 159 | 31 | 24 | 7 | 0.281 |
| C2 | `angle` | `z_noise5` | `z_noise0` | 91 | 100 | 20 | 29 | −9 | 0.852 |
| C3 | `angle` | `z_noise0` | `oracle_noise0` | 100 | 159 | 23 | 82 | −59 | 1.000 |
| C4 | `angle` | `z_noise5_corrected` | `z_noise0_corrected` | 84 | 88 | 22 | 26 | −4 | 0.702 |
| C1 | `shape` | `oracle_noise5` | `oracle_noise0` | 120 | 114 | 43 | 37 | 6 | 0.334 |
| C2 | `shape` | `z_noise5` | `z_noise0` | 123 | 118 | 40 | 35 | 5 | 0.363 |
| C3 | `shape` | `z_noise0` | `oracle_noise0` | 118 | 114 | 70 | 66 | 4 | 0.422 |
| C4 | `shape` | `z_noise5_corrected` | `z_noise0_corrected` | 93 | 107 | 29 | 43 | −14 | 0.904 |

Per reading (descriptive; the same file, scopes `<run>/<control>`), C2's `delta` discordance concentrates at the
translation control:

| reading | B only | A only | p |
|---|---|---|---|
| δ = 0.1, translation | 13 | 0 | 0.0001 |
| δ = 0.05, translation | 8 | 0 | 0.004 |
| δ = 0.1, anchor | 9 | 9 | 0.593 |
| δ = 0.05, anchor | 11 | 9 | 0.412 |

The same holds, more weakly, for the corrected C4: 9 / 3 and 8 / 1 at translation, 9 / 9 and 6 / 8 at the anchor.

## 4. Rank tables

Anchor units by selected rank, with rejections `delta` / `angle` / `shape` (`diagnosis_rank_table.csv`):

| arm | run | rank 1 | rank 2 | rank 3 | rank ≥ 4 |
|---|---|---|---|---|---|
| `z_noise5` | δ = 0.1 | — | 366: 25 / 15 / 19 | 128: 9 / 9 / 9 | 6: 0 / 0 / 0 |
| `z_noise5` | δ = 0.05 | — | 407: 28 / 20 / 21 | 92: 4 / 4 / 6 | 1: 0 / 0 / 0 |
| `latent_stage_cv` | δ = 0.1 | — | 344: 22 / 35 / 27 | 65: 0 / 6 / 10 | 91: 8 / 4 / 9 |
| `latent_stage_cv` | δ = 0.05 | — | 383: 24 / 34 / 27 | 56: 2 / 6 / 7 | 61: 6 / 2 / 5 |
| `latent_stage_floor1` | δ = 0.1 | 56: 6 / 0 / 2 | 297: 18 / 32 / 22 | 60: 0 / 5 / 9 | 87: 7 / 4 / 8 |
| `latent_stage_floor1` | δ = 0.05 | 60: 5 / 0 / 1 | 330: 22 / 32 / 22 | 54: 2 / 6 / 7 | 56: 4 / 2 / 4 |
| `latent_stage_rank2` | δ = 0.1 | — | 500: 30 / 42 / 34 | — | — |
| `latent_stage_rank2` | δ = 0.05 | — | 500: 31 / 41 / 32 | — | — |

- **`latent_stage_cv` by rank.**
  - Its `angle` excess sits at rank 2: 35 / 344 = 0.102 and 34 / 383 = 0.089. The tail rejects less (4 / 91 and
    2 / 61 at rank ≥ 4).
  - Its `shape` excess sits in the tail: 10 / 65 and 7 / 56 at rank 3, 9 / 91 at rank ≥ 4, against 27 / 344 and
    27 / 383 at rank 2.
- **Forcing rank 2** therefore removes the `shape` excess (34 / 32 at the anchor) but not the `angle` excess (42 / 41).
- **Floor 1.** It selects rank 1 in 56 / 500 (0.112) and 60 / 500 (0.120) anchor units (`"floor1"` in the decision
  file). Those units never reject `angle`, the two-point null of the [rank-floor check](rank-floor-2026-10-23.md).
- **`z_noise5`** is at rank 2–3 in 494 and 499 anchor units. Its `delta` rejections are spread over ranks 2 and 3 in
  proportion (0.068 and 0.070 at δ = 0.1), so its `delta` excess is not a rank effect.

## 5. What the arms show

**The case study's `delta` excess is residual-df shrinkage.** The mechanism the proposal named predicts three things,
and the run shows all three:
- **The uncorrected rate rises with the reduced model's column count in both families.**
  - With `z`: 0.0545, 0.0595, 0.066 at p_red 5, 7, 10.
  - Under the oracle: 0.0525, 0.0545, 0.060 at p_red 4, 6, 9.
  - The proposal's normal-tail estimate was 0.067 at six columns and 0.056 at none.
- **The corrected rate is flat and at or below α.**
  - With `z`: 0.049, 0.048, 0.0545.
  - Under the oracle: 0.048, 0.0475, 0.0505.
  - No corrected arm exceeds the bound on `delta` in any reading; `z_noise5_corrected`'s largest is 30 / 500.
- **`angle` is not inflated.** The `z` arms read 0.042–0.050.

C2 is significant at 0.015. C1 points the same way but stops at 0.065 (net 15 against C2's 23). C4 keeps a residual
net of 11 (p = 0.156), so the uniform κ removes most, though not visibly all, of the column trend. C3 is flat (p = 0.39,
net 4), so estimating a covariate confounded with the group × stage cells costs no `delta` size here.

The C2 discordance comes from the translation control: 21 readings against 0, while the anchor readings balance at
20 against 18. That fits a null narrowed by the reduced model's columns, which a constant location offset in group B
leaves otherwise untouched. The parent's single excess (translation `delta` at δ = 0.1, 36 / 500) belongs to this
pattern.

**The correction costs some `angle` and `shape` rejections.** `z_noise5_corrected` reads `angle` 0.042 and `shape`
0.0465, against 0.0455 and 0.0615 uncorrected. A uniformly widened null is conservative for the two statistics that
were not inflated. That is a power question for the follow-up change (design D5 deferred the leverage-adjusted
alternative to it).

**The oracle is not a size-α reference for `angle`.**
- Removing `z`'s injection exactly from the latent scores, while leaving `z` out of the model, gives `angle`
  0.0725–0.083 in every reading, corrected or not. Every `z` arm holds `angle`. C3 on `angle` is 23 against 82 in the
  opposite direction.
- The basis is fitted on features that contain `z`, and only a `z` column in both models absorbs whatever that basis
  selected along `z`. That is a candidate explanation, not a tested one: the smoke probes (`design.md`, "Smoke and
  shard plan") found no `z`-aligned content in `Y0`'s within-cell or interaction pattern on 10–20 units. The
  mechanism stays open.
- The smoke had flagged the oracle rows for `delta` (5–7 / 20). At R = 500 the oracle's `delta` is 0.0525–0.060, so
  that was noise at R = 20. The `angle` inflation is real.
- The `delta` verdict reads only paired contrasts, and C1 is within the oracle family. Neither depends on the oracle
  being size-α. Only C3's absolute meaning does, and C3 is flat on `delta`.

**The `stage` latent adjustment fails on `angle` at rank 2, not through its rank tail.**
- With the rank forced to 2 and no cross-validation, which removes both the rank tail and the per-training-set
  nuisance coefficients, anchor `angle` still reads 0.084 / 0.082.
- Under CV, its rank-2 units reject `angle` at 0.102 / 0.089. The parent's `latent_stage` `angle` excess is
  therefore a property of the residualized basis.
- `shape` behaves differently: its excess came with the CV rank's tail, and at rank 2 it holds (0.068 / 0.064).
- Floor 1 changes little: anchor `angle` 0.082 / 0.080 and `shape` 0.082 / 0.068, with 11–12 % of units at rank 1.

## 6. Verdicts and consequences

Read verbatim from `case_study_diagnosis_decision.json`.

**Gate.** `"reproduction": {"status": "holds"}`, with 3,000 / 0 compared / mismatching per arm and run. No run was a
smoke run, so `"indicative": false`.

**Covariate-adjusted `delta`: `residual_df_shrinkage`.**
- Conditions:
  - `"column_contrast_significant": {"C1": false, "C2": true}`
  - `"corrected_case_study_model_holds": true`
  - `"corrected_column_contrast_significant": false`
  - `"collinear_contrast_significant": false`
- p-values C1 0.0645, C2 0.0153, C3 0.3863, C4 0.1561.
- **For MOTCO.** A change offering the residual-df corrected permutation:
  - opt-in, with every output byte-identical when it is absent;
  - its own Type I and power check, comparing the uniform κ with the leverage-adjusted residuals of design D5, and
    reporting the `angle` and `shape` cost seen here.
- **For `sead_sexdiff`.**
  - The primary stays withheld until that change is released and validated; then the primary is re-issued under it.
  - Until then, the primary's `delta` p-values read against a realized size of about 0.066 at six covariate columns.
    That is now attributed to the reduced model's column count, not to the covariate's confounding.
  - Its `angle` and `shape` p-values are size-controlled (0.042–0.068 in `z_noise5`).
  - Nothing in `sead_sexdiff` changes under this report.

**`stage` latent adjustment: `residualized_basis`.**
- `"rank2_not_holding"` lists anchor `angle` at δ = 0.1 (0.084) and δ = 0.05 (0.082). `latent_stage_rank2` holds every
  other reading.
- **For `sead_sexdiff`.** The `stage` sensitivity arm stays inadmissible at this configuration, at a fixed rank of 2
  too. Its `angle` rejections are not interpretable there. A rank-2 `stage` arm would hold `shape`, but the verdict
  is on the arm, so it is not adopted.
- **For MOTCO.** Separating the per-training-set coefficients from the residualized basis would need a core switch
  (design, Risks). It is not proposed, because the forced-rank arm already has neither the rank tail nor CV.

## 7. Limitations and scope

The parent check's four limitations carry over unchanged:

- **Astrocyte primary cohort sizes only.** `[[11, 10, 28], [9, 10, 12]]` with n = 80. The microglia cohort (n = 76)
  and the other cohorts' tables were not run. The shrinkage scales with p_red / n, so a cohort with fewer donors at
  the same covariate set is expected to be more liberal, not less.
- **Additive-linear covariate effects.** `z` enters the features additively along one fixed direction and the
  trajectory model linearly.
- **One confounded covariate beside five independent ones.** The five noise covariates are i.i.d. normals,
  uncorrelated with one another and with the design, unlike the case study's age, education, PMI, RIN and library
  sizes. The column-count mechanism does not depend on their joint distribution, but its size at the real covariates
  is inferred.
- **Structural-analogue blocks.** These are InterSIM-calibrated methylation and expression analogues (498 measured
  features), not the case study's TMM log2-CPM blocks of 5,000 features each.

And those of the diagnosis itself:

- **`latent_stage_rank2` confounds two things.** It removes the CV rank's tail and the per-training-set nuisance
  coefficients together. Under `residualized_basis` neither is the cause of the `angle` excess, but the arm cannot say
  which of them carried the parent's `shape` excess. The rank table places that excess in the tail.
- **The oracle's absolute rates.** The oracle arms are a valid within-family reference for the column contrast C1, not
  a size-α baseline for `angle` (§5).
- **No new replicates.** The diagnosis reuses the parent's 500 replicates per strength, so the parent's own excess
  (two rejections over the bound) is explained rather than re-measured. The R = 1,000 replication was not needed under
  this verdict and stays deferred.
- **The correction is minimal.** The uniform κ restores the unbiased residual variance but not the leverage structure,
  and C4's residual net of 11 (p = 0.156) leaves room for that. The follow-up change measures both corrections' size
  and power.
- **No package change.** Nothing under `src/motco/` changed: the correction is computed in the script from the public
  `RRPP`, and MOTCO stays at 0.12.0. The case study's pin is untouched.
- **What stays on `ing`.** The 45,000-row `diagnosis_pvalues.csv` (2.5 MB per run), `merged.jsonl` and `records/` stay
  there, named with their sha256 in each `PROVENANCE.txt`. `case_study_diagnosis.py merge` regenerates every committed
  table from `records/`.
