# Predictive-subspace size excess — diagnosis report (2026-11-13)

**Change:** `predictive-size-excess-diagnosis` (no package change; MOTCO 0.24.0). **Reads:**
`results/predictive-size-excess-2026-11-13/` (`predictive_size_excess_decision.json`; the rate, paired and cell
tables; `accounting.json`, `reproduction.json`, `preflight.json`, `settings.json`; `PROVENANCE.txt`, `sacct.txt`).
**Parent:** the [predictive-subspace pilot](predictive-subspace-2026-11-12.md). **Concept page:**
[measurement visibility](../measurement-visibility.md#the-measured-dimension-fit-rank-against-predictive-subspace).

The predictive-subspace pilot read `keep_production` on one finding. On the SEA-AD plasmode anchor (500
permuted-label replicates of the astrocyte export) all three predictive coordinate systems rejected `delta` in 65–68
units and `shape` in 38–39, where production rejected 21 and 13. No simulated design reproduced the excess.

Roadmap open item 2 named a guess. The two predictive directions of the rank-5 fit hold its in-sample stage fit
(stage R² 0.84 against a held-out AUROC of 0.636), and the reduced model's residuals in them are not exchangeable
across stages, which RRPP's draws assume. Two checks test it on the pilot's own units, labels and RRPP seeds, under a
decision fixed before any unit ran (`design.md` D8):

- **Check W** keeps the coordinates and changes the null. RRPP's draws are restricted within stage: the reduced
  model's residuals are permuted among the donors of each stage, never across stages.
- **Check O** keeps the null and changes the coordinates. The fitted stage indicators are computed out of fold:
  each donor's row comes from a stage-PLS model fitted without it.

**In one paragraph.** **The verdict is `both_restore`.** Each repair alone brings the predictive coordinates' anchor
rates within the size rule the pilot applied, so the excess needs both the unrestricted draws and coordinates fitted
on the donors they measure.

- **Within-stage draws remove the `delta` excess entirely.** 37–38 units are rejected only under the unrestricted
  draws, and none only under the within-stage ones; `delta` falls from 0.130–0.136 to 0.056–0.062. They remove
  little of the `shape` excess. `shape` falls from 0.076–0.078 to 0.064–0.066 (5–6 discordant units, all one way),
  so it holds the rule by one or two counts and stays at about 2.5 times production's 0.026.
- **The out-of-fold coordinates read 0.016 / 0.014 / 0.034**, but at the anchor they carry almost no stage signal:
  the out-of-fold fit explains a median 0.022 of the stage indicators' variance, against 0.84 in sample.
- **On the real-noise plant** the within-stage arms keep `delta` above production's at e = 0.5 (0.64–0.66 against
  0.51), but production gains as much under the same draws (0.63).
- **Without covariates** the predictive excess stays, so the six covariates and their `df` correction do not cause
  it.
- **The rule's consequence.** A follow-up change may propose within-stage draws as an opt-in option and read their
  size and power before `predictive_orthonormal` is re-read. The pilot's `keep_production` stands.

## 1. Configuration and provenance

- **The units.** The pilot's 800 plasmode units (`predictive_subspace_check.plasmode_settings` / `plasmode_units`,
  imported unchanged): the anchor at the `paper` stream's replicates 0–499, and the production selection's real-noise
  span plant at e ∈ {0.25, 0.5, 1.0} at replicates 0–99. All run at 999 permutations with rejection at p ≤ 0.05.
  Each unit carries the pilot's label draw (within-stage, base seed 2034), plant and RRPP seed. The configuration is
  the case study's: six covariates, the `df` residual correction, net displacement and rank floor 2. The export is
  the SEA-AD MTG astrocyte `cohort_30` (RNA and ATAC, 5,000 features each, 80 donors, cells F 11 / 10 / 28 and
  M 9 / 10 / 12), read in place through the plasmode's contract.
- **One fit per unit.**
  - An anchor unit fits at the plasmode's committed fidelity rank of 5 without cross-validation. The anchor does not
    move the label-free features, so every anchor unit's CV selects that rank.
  - A plant unit runs the production CV once, and every other arm is forced at the rank it selects (4–14; medians 7, 7
    and 5 at e = 0.25, 0.5 and 1.0).
  - The pre-flight required the cross-validated fit on anchor replicates 0 and 1 to select rank 5 at AUROC 0.635556,
    and to give the forced call's statistics and p-values exactly. It passed on the workstation and on `ing`.
- **The fourteen arms** (`design.md` D7), each at the unit's RRPP seed with serial draws:

  | arm | coordinates | draws | covariates | role |
  |---|---|---|---|---|
  | `pls_cv` | scores | unrestricted | six | reference |
  | `predictive`, `predictive_orthonormal`, `predictive_fitted` | in-sample predictive | unrestricted | six | reference |
  | `predictive__within_stage`, `predictive_orthonormal__within_stage`, `predictive_fitted__within_stage` | in-sample predictive | within stage | six | check W, decisive |
  | `fitted_out_of_fold` | out-of-fold fitted indicators | unrestricted | six | check O, decisive |
  | `pls_cv__within_stage` | scores | within stage | six | descriptive |
  | `fitted_out_of_fold__within_stage` | out-of-fold fitted indicators | within stage | six | descriptive |
  | `<reference arm>__no_covariates` (four) | as named | unrestricted | none | descriptive |

  - **Within-stage draws.** Per draw, each stage's residual rows are permuted among themselves, in declared stage
    order, from `default_rng(the unit's RRPP seed)`. Everything else is `RRPP`'s: the reduced fit, the `df` scale and
    `estimate_difference` per draw. The driver's loop equals `RRPP` exactly when every row is in one stratum (a
    test). One restricted sequence per unit is shared by every within-stage arm.
  - **Out-of-fold coordinates.** Four stage-stratified, shuffled folds per replicate, seeded from the driver's name,
    base seed 2037, the replicate and `folds`. For each fold, `fit_plsda_model` is fitted on the other three folds'
    rows of the analysis's PLS input features, at the unit's rank and scaling. The held-out rows then get
    `transform(X_held) · y_loadingsᵀ`, the fitted stage indicators in that model's response units. With every donor in
    training, the same computation returns `predictive_fitted` to 1e-8 (a test). Each training set holds 15 / 15 / 30
    donors by stage.
- **Checks on every unit.**
  - The pilot's five consistency checks.
  - `same_fit_covariate_free`: each covariate-free arm fitted the production model.
  - `path_identity`: the driver's inference path, under unrestricted draws, returns each reference arm's statistics
    to 1e-10 and its p-values exactly.
- **Code and host.** Revision `7046672` on `proposal/predictive-size-excess-diagnosis`, run from a clean detached
  worktree on `ing` (partition `512x1024`, nodes n1 and n3). Python 3.11.16, numpy 2.3.5, scipy 1.16.3,
  scikit-learn 1.8.0 and pandas 2.3.3, with BLAS pinned to one thread. Preparation job 938594, array 938598
  (100 shards) and finishing job 938732; every allocation completed with exit 0:0. 10.92 recorded core-hours, against
  a projection of about 11.5 and a cap of 20 (`PROVENANCE.txt`).

## 2. Accounting and reproduction

- **Units.** 800 of 800 complete: no failed unit, no failed arm, no inconsistent unit and no duplicate record.
- **Checks.** Every applicable check held on every unit. The largest deviations were `path_identity` 0, `same_fit` 0,
  `same_fit_covariate_free` 0, `rank` 0 and `fitted_labels` 5.3e-14. Every covariate-free arm's latent matrix equals
  its reference arm's exactly.
- **Reproduction: passed.** Each of the four reference arms' available and rejected counts, in each of the four
  cells, for `delta`, `angle` and `shape`, equals the pilot's committed `rate_table.csv`: 48 of 48 rows. The forced
  anchor fit therefore reproduces the pilot's cross-validated anchor counts exactly (21 / 24 / 13 for `pls_cv`, 68 / 11
  / 38 for `predictive`, 67 / 9 / 38 for `predictive_orthonormal`, 65 / 10 / 39 for `predictive_fitted`).

## 3. The decision

Rule version 1 grades each decisive arm on its three anchor readings:

- **`restored`:** the arm holds the stage-span paper grade's count-plus-ceiling rule (`judge_rates`, m = 3). This is
  the rule and the family the pilot's S1 failed. At R = 500 a reading exceeds the bound α + 2·SE = 0.0695 above 34
  rejections, the ceiling is 36, and one exceedance is allowed.
- **`reduced`:** not restored, and the one-sided exact paired p-value for fewer rejections than the arm's
  counterpart is at most 0.025 on `delta` or `shape`.
- **`unchanged`:** otherwise.

Check W takes the weakest grade of its three arms, and check O the grade of its one arm.

| arm | `delta` | `angle` | `shape` | exceedances | grade | pilot's counterpart |
|---|---|---|---|---|---|---|
| `predictive__within_stage` | 31 (0.062) | 17 (0.034) | 33 (0.066) | 0 | `restored` | 68 / 11 / 38 |
| `predictive_orthonormal__within_stage` | 29 (0.058) | 17 (0.034) | 32 (0.064) | 0 | `restored` | 67 / 9 / 38 |
| `predictive_fitted__within_stage` | 28 (0.056) | 17 (0.034) | 33 (0.066) | 0 | `restored` | 65 / 10 / 39 |
| `fitted_out_of_fold` | 8 (0.016) | 7 (0.014) | 17 (0.034) | 0 | `restored` | 65 / 10 / 39 (`predictive_fitted`) |

Check W is `restored`, and so is check O. **The verdict is `both_restore`.** The gate passed: the reproduction
passed, no unit was inconsistent and the pre-flight passed.

The consequence fixed with the rule: "The unrestricted draws are what fails in the predictive subspace on this
export, and the out-of-fold reading confirms that the in-sample fit is where the excess lives. A follow-up change
may propose within-stage draws as an opt-in option and must read their size and power on the pilot's simulated
designs, the weak design and this anchor before `predictive_orthonormal` is re-read."

## 4. Check W: what the within-stage draws repair

Paired against each arm's counterpart (the same coordinates under unrestricted draws) on the 500 anchor units. b
counts the units only the counterpart rejects, and c the units only the within-stage arm rejects:

| arm | `delta` b / c | `delta` p (fewer) | `shape` b / c | `shape` p (fewer) |
|---|---|---|---|---|
| `predictive__within_stage` | 37 / 0 | 7.3e-12 | 5 / 0 | 0.031 |
| `predictive_orthonormal__within_stage` | 38 / 0 | 3.6e-12 | 6 / 0 | 0.016 |
| `predictive_fitted__within_stage` | 37 / 0 | 7.3e-12 | 6 / 0 | 0.016 |

- **`delta`.** The within-stage null removes the whole excess, and only that. Every unit the restricted draws stop
  rejecting was one the unrestricted draws rejected, and none goes the other way. Against production's unrestricted
  rates the within-stage arms still read slightly more (17–19 units against 9–10; two-sided p
  0.09–0.25).
- **`shape`.** It barely moves: 5–6 units, all one way. The arms hold the rule at 32–33 rejections against the
  exceedance threshold of 34. Against production `shape` stays higher: 26 units against 6–7, two-sided p 0.0005–0.0013.
  Production itself reads 19 under the same within-stage draws (§7).
- **`angle`** rises from 9–11 to 17 and stays well inside the bound.

So `restored` is firm for `delta` and marginal for `shape`. The rule grades both alike; the follow-up has to read
`shape` again.

## 5. Check O: what the out-of-fold coordinates show

`fitted_out_of_fold` reads 8 / 7 / 17 on the anchor. Against `predictive_fitted` that is 64 / 7 discordant units on
`delta` and 38 / 16 on `shape`. Against production's unrestricted rates it rejects less `delta` (20 / 7) and less
`angle` (21 / 4), with `shape` close (13 / 17).

The reading is weaker evidence than its grade suggests. At the anchor the out-of-fold coordinates carry almost no
stage signal. The out-of-fold fit explains a median 0.022 of the stage indicators' variance, against 0.84 in sample
on every anchor unit. Its `angle` null has a median 95th percentile of 94.9°, so the net displacement it measures is
close to noise.

What the arm shows is that the excess needs coordinates fitted on the donors they measure. It does not show that a
cross-fitted measurement is usable at this design. Under within-stage draws the same coordinates read 23 / 18 / 23
(descriptive), also inside the rule.

## 6. The plant cells (descriptive)

`delta` rejection rates at R = 100. The paired tests are one-sided, for a gain against `pls_cv`:

| arm | e = 0.25 | e = 0.5 | e = 1.0 | e = 0.5 against `pls_cv` (b / c, p) |
|---|---|---|---|---|
| `pls_cv` | 0.24 | 0.51 | 1.00 | — |
| `predictive`, `predictive_orthonormal`, `predictive_fitted` | 0.40–0.42 | 0.79–0.82 | 1.00 | 1–2 / 30–32, ≤ 1.2e-7 |
| `predictive*__within_stage` | 0.28–0.29 | 0.64–0.66 | 1.00 | 8–9 / 22–23, 0.005–0.015 |
| `pls_cv__within_stage` | 0.24 | 0.63 | 1.00 | 0 / 12, 0.00024 |
| `fitted_out_of_fold` | 0.81 | 0.99 | 1.00 | 0 / 48, 3.6e-15 |
| `fitted_out_of_fold__within_stage` | 0.82 | 0.99 | 1.00 | — |

- **The within-stage arms keep part of the real-noise gain.** At e = 0.5 they read 0.64–0.66, above production's
  0.51 (the decision's statement for each `restored` arm), and below their unrestricted counterparts' 0.79–0.82
  (14–16 units, all one way).
- **Production gains under the same draws.** It reads 0.63 at e = 0.5 (12 units, all one way). The within-stage
  draws keep the pooled stage means at their fitted values, which narrows the `delta` null against a magnitude edit
  of one group's configuration. Within-stage draws against within-stage draws, the predictive coordinates' advantage
  is therefore small: 0.64–0.66 against 0.63 at e = 0.5, and 0.28–0.29 against 0.24 at e = 0.25. No paired test of
  that pair was predeclared, so these are rates only.
- **The out-of-fold arm reads 0.81 and 0.99, above every in-sample arm.** The design expected the opposite, a lower
  bound. The plant is defined by the full sample's stage means. A plausible reading, not tested here: each fold
  model is fitted on 60 of the same 80 donors, so it recovers the planted direction, while a held-out donor's own
  noise is not fitted into its coordinates. On that reading the out-of-fold plant rate reflects how the plant is
  built, not real-noise power. The arm also leaks `shape` there (0.16 at e = 0.5).

## 7. Production under within-stage draws (descriptive)

`pls_cv__within_stage` reads 19 / 27 / 19 on the anchor, inside the rule. Against `pls_cv` the paired counts are 4 / 2
on `delta`, 0 / 3 on `angle` and 0 / 6 on `shape` (two-sided p 0.69, 0.25 and 0.031). Production holds size under
either null on this anchor. On the plant it gains `delta` power under the within-stage draws (§6) and reads more
off-target `shape` (0.11 at e = 0.25 and e = 0.5, against 0.09 and 0.04).

## 8. The covariate-free arms (descriptive)

| arm | `delta` | `angle` | `shape` | against the six-covariate counterpart, `delta` (b / c, two-sided p) |
|---|---|---|---|---|
| `pls_cv__no_covariates` | 41 (0.082) | 14 | 16 | 3 / 23, 0.0001 |
| `predictive__no_covariates` | 63 (0.126) | 14 | 42 | 12 / 7, 0.36 |
| `predictive_orthonormal__no_covariates` | 61 (0.122) | 14 | 41 | 13 / 7, 0.26 |
| `predictive_fitted__no_covariates` | 60 (0.120) | 15 | 44 | 14 / 9, 0.40 |

- **The predictive excess does not come from the covariates.** Without them the three arms read the same `delta` and
  `shape` excess as with them: no paired difference has p below 0.26.
- **Production without covariates exceeds on `delta`.** It reads 0.082, above the bound of 0.0695, with 23 units
  against 3. The `df` correction then follows the four-column reduced model. On this anchor production's `delta`
  size depends on having the covariates in the trajectory model. The case study runs with them; this reading is a
  lead about the covariate-free configuration, not about the case study.

## 9. Dispersion by stage and the stage fit (descriptive)

On each unit, each stage's mean squared residual norm under the six-covariate reduced model, divided by the pooled
value. The table gives medians over the 500 anchor units; stages are 0 / 1 / 2 in declared order, with 20 / 20 / 40
donors.

| coordinates | stage 0 | stage 1 | stage 2 | max / min |
|---|---|---|---|---|
| scores | 0.805 | 0.921 | 1.137 | 1.41 |
| `predictive` | 1.121 | 1.178 | 0.852 | 1.39 |
| `predictive_orthonormal` | 1.120 | 1.278 | 0.802 | 1.59 |
| `predictive_fitted` | 1.091 | 1.252 | 0.830 | 1.51 |
| out of fold | 0.841 | 0.564 | 1.293 | 2.32 |

- **The size of the inequality does not separate the predictive coordinates from the scores** (1.39–1.59 against
  1.41).
- **Its direction does.** In the scores the largest stage is the widest. In the predictive coordinates it is the
  narrowest (0.80–0.85), and the two 20-donor stages are wider (1.09–1.28). The in-sample fit tightens the stage with
  the most donors along the directions it fitted.
- **What the unrestricted draws do with it.** They move residuals between stages of unequal dispersion, and the
  within-stage draws do not. The measured `delta` excess is consistent with that. The table describes the pattern; it
  does not test the mechanism.
- **The stage fit.** In-sample stage R² is 0.840 on every anchor unit, against a median out-of-fold share of 0.022.
  On the plant cells the two readings are 0.95 / 0.27, 0.95 / 0.40 and 0.87 / 0.49 at e = 0.25, 0.5 and 1.0. The
  plant adds stage signal that survives out of fold.

## 10. Null widths (descriptive)

Median `angle` null 95th percentile on the anchor:

| arm | degrees |
|---|---|
| `pls_cv` | 40.2 |
| `pls_cv__within_stage` | 38.8 |
| `predictive`, `predictive_orthonormal`, `predictive_fitted` | 11.0–12.1 |
| their within-stage arms | 9.3–10.2 |
| `fitted_out_of_fold` | 94.9 |
| `fitted_out_of_fold__within_stage` | 74.1 |

The predictive coordinates' `angle` null is narrow, as in the pilot. Within-stage draws narrow it further, and
their `angle` rejections rise from 9–11 to 17, still well inside the bound.

## 11. What follows

- **The verdict's consequence.** A separate change may propose within-stage draws as an opt-in option of `RRPP`, the
  harness and `motco analyze`, byte-identical when unset, and must read their size and power on the pilot's
  simulated designs, the weak design and this anchor. Only its verdict could reopen `predictive_orthonormal`, under
  the admissibility rule the roadmap requires.
- **What that change should carry from this reading.**
  - **`shape` stays near the bound** under within-stage draws on this anchor (0.064–0.066), so its size reading
    cannot be assumed.
  - **Production gains as well.** Production gains plant power under within-stage draws (0.51 → 0.63 at e = 0.5), so
    any comparison of measurement spaces must pair them under the same draws. The predictive coordinates' advantage
    over production then shrinks to a few points on the plant.
  - **The null is a design question, not only a predictive-subspace one.** Within-stage draws assume exchangeability
    within stage only, which is weaker than `RRPP`'s assumption, and they keep the pooled stage means at their fitted
    values.
- **The pilot's `keep_production` stands.** The predictive values stay opt-in diagnostics. No verdict here validates
  a test or changes a default.
- **Open item 2 keeps its prediction for ADNI.** The guess was confirmed on this export, so the predicted excess
  should shrink as donors per feature grow. Reading it needs ADNI's export contract.
- **A covariate-free lead.** Production without covariates exceeds on `delta` at this anchor (§8). Whether that holds
  beyond one fit is unread.

## 12. Limitations

- **One export and one fit.** The 500 anchor units share one latent matrix and differ only in the group labels, so
  every anchor rate is conditional on that fit (as the pilot's were). The out-of-fold arm averages over fold splits,
  not over donors.
- **The within-stage null mirrors how the plasmode draws its labels.** The plasmode permutes the labels within stage,
  so the restricted draws match the anchor's null by construction. Whether they hold size where residuals are
  exchangeable across stages, and what they cost in power, is the follow-up's simulation reading.
- **`shape` holds by one or two counts.** At R = 500 the within-stage arms sit at 32–33 against a threshold of 34. A
  rerun on new label draws could put one of them over it.
- **The out-of-fold arm's `restored` comes with almost no stage signal** (§5), and its plant reading reflects the
  plant's construction (§6).
- **The plant cells are R = 100 and descriptive.** No paired test between the within-stage predictive arms and
  production under the same draws was predeclared.
- **The rank is chosen on all donors.** It is one integer, selected by held-out AUROC and pinned at 5 on the anchor.
