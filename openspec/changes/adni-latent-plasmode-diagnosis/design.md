# Design

## Context

See `proposal.md` (Why) and `notes.md` (the exploratory pass). The parent run's units each record, per space and reading, a numerator `N` (the halves' cross-product inside the space) and a denominator `D` (in total), each averaged over 20 splits. A calibration unit also records the truths `TN` and `TD`: the planted configuration's energy inside the space and in total, in the unit's coordinates. A cell holds 24 synthetic cohorts × 25 units (calibration) or 500 units of the one cohort (ADNI). The parent's estimate is `ΣN / ΣD` over a cell's units, per cohort in the calibration. Its gate judged the mean over cohorts of that ratio's error.

## Goals / Non-Goals

**Goals**

- Say whether the gate's bias lies in the estimator's parts or in the ratio.
- Say whether the spread would shrink with more units or splits of one cohort.
- Say whether any interval for ADNI's capture is bounded at the run's points.

**Non-Goals**

- A new unit, a new calibration strength or a larger synthetic cohort.
- Changing the parent's verdict. `gate_failed` stands; the diagnosis says what follows from it.
- A rule version 2 for the capture. If the reading allows one, it is its own change.

## Decisions

### D1. No new units; the records are bound

- **Inputs.** `<records-dir>/<profile>/merged.jsonl` for both profiles, read where the run wrote them.
- **Bindings.** The parent's five committed files (both settings, the calibration table, the cohort's capture table, the decision) are bound by sha256 in `settings.json` and re-checked by `run` and `decide`. Every record must carry its profile's committed settings hash, and each expected unit must be present exactly once. Each records file's sha256 is written to `reproduction.json`.
- **Reproduction.** The parent's `calibration_rows` and `capture_rows`, written by its CSV writer, must equal the committed files byte for byte. If either differs, the verdict is `reproduction_failed`.

### D2. F1, the parts

For each calibration cell with a truth (strengths 0.67 and 1) and each space and reading:

- the mean of `N − TN` and of `D − TD` over the cell's units;
- a standard error from the 24 cohort means of those differences (units of one cohort are dependent);
- `z` = mean / standard error.

### D3. F2, the aggregate

The ratio of sums over all of a cell's 600 units against the truth's ratio of sums. Its standard error and 95 % interval come from 2,000 cluster-bootstrap resamples of cohorts (base seed 2042, one stream per cell). If the parts are unbiased and the aggregate is on its truth, the per-cohort bias is the nonlinearity of a ratio with a noisy denominator.

### D4. F3, where the spread lives

- `W` is the within-cohort covariance of `(N, D)`, pooled over cohorts with ddof 1 per cohort.
- `C` is the covariance of the 24 cohort means.
- `B = C − W / m`, with `m` = 25 units per cohort.
- The between share of `var(D̄)` is `B_DD / C_DD`. When it is near 1, a cohort's mean does not move with more units or splits.
- Recorded beside it:
    - the per-cohort ratio's SD and its robust SD (IQR / 1.349);
    - the cohorts whose mean `D` is below a quarter of the cell's;
    - the delta-method bias and SD at `m`, at 500 units and as units grow without bound. The delta method is a guide only: it fails where the denominator's coefficient of variation is large, which the table also records.

### D5. F4, ADNI's own interval

- **ADNI's inputs.** For each point, the mean `(N̄, D̄)` over the cohort's 500 units, and their covariance `W_ADNI`.
- **The covariance of the mean.** It is `V = f·B + W_ADNI / 500`. ADNI's between-cohort covariance cannot be observed from one cohort, so `B` is taken from the calibration at the same point, space and reading. There are four stand-ins:
    - `B` at strength 0 or at strength 0.67;
    - each taken `raw` (`f = 1`) or `scaled` (`f` is ADNI's within-cohort variance of `D` over that strength's).
- **The primary stand-in** is strength 0.67, raw, the strength the parent's interval used.
- **The interval.** Fieller's 95 % interval for `N / D` from `(N̄, D̄, V)` is `bounded` only when `t = D̄ / √V_DD` exceeds 1.96. Otherwise it is the line outside two roots (`exclusive`) or the whole line.
- **Why Fieller.** It is the interval for a ratio that stays valid when the denominator's uncertainty is not small. A bounded Fieller interval is necessary for any rule that states a capture with an interval.

### D6. The reading (rule version 1)

Written after the exploratory pass (`notes.md`). It is therefore a reading with its thresholds stated, not a predeclared test, and it does not amend the parent's rule.

- **F1** holds when `|z|` of both parts is below 3 at every gate cell (`pls_cv`, total).
- **F2** holds when `|z|` of the aggregate's bias is below 3 at every gate cell.
- **F3** holds when the between share of `var(D̄)` at strength 0.67 is at least 0.5 at every point.
- **F4** holds when ADNI's interval is unbounded at every point under the primary stand-in. Agreement over the four stand-ins is recorded.

| verdict | condition | what it means |
|---|---|---|
| `reproduction_failed` | a table differs | nothing is read |
| `estimator_biased` | F1 or F2 fails | the estimator needs repair before any rule version 2 |
| `resolvable` | F4 fails | a rule version 2 with Fieller's interval could read where it is bounded |
| `units_limited` | F3 fails | more units per point could narrow the spread |
| `information_limited` | otherwise | the estimator is sound, its spread is the cohort's, and ADNI's program does not clear it: no rule version 2 on ADNI |

The limit 3 is set for 16 tests (8 cells × 2 parts) in F1 and 8 in F2. The limit 0.5 separates "mostly between" from "mostly within".

### D7. Outputs

- `component_table.csv` (F1, F2), `spread_table.csv` (F3 and the per-cohort spread) and `adni_table.csv` (F4): one row per strength or stand-in, point, space and reading.
- `reproduction.json` and `settings.json`.
- From `decide`, `adni_latent_plasmode_diagnosis.json`: the verdict, the findings per cell, the rule, the parent's verdict and the inputs' hashes.
- No per-unit and no per-cohort value is written.

## Risks / Trade-offs

- **[The reading was written after seeing the numbers]** → It is labelled a reading in the rule, the decision file and the report. Its thresholds are stated with their reasons. The parent's verdict, the only predeclared one, is unchanged.
- **[ADNI's between-cohort covariance is a stand-in]** → Four stand-ins are recorded, and F4 is judged on the one the parent used. The report gives `t` under each.
- **[24 cohorts estimate a variance roughly]** → `B` carries about 30 % relative error, and the report states this. F4's margin is read beside it.
- **[The delta method fails near a zero denominator]** → Only F2 (bootstrap) and F4 (Fieller) enter the reading. Delta-method values are recorded for description.
