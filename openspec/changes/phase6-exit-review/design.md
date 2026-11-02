# Design — phase6-exit-review

## Context

See `proposal.md` (Why) for the motivation. This section gives the state the review decides on.

- **The gate and the Work list** (`docs/roadmap.md`, Phase 6). Phase 6's purpose is "demonstrate that the validated method answers a meaningful biological question end to end". Its exit gate has two parts:
  - **G1**: trace each global trajectory claim to stable feature/module drivers.
  - **G2**: reproduce the result from documented inputs and configuration.

  The Work list has six bullets. One of them is conditional: "For significant orientation, produce observed and PLS-reconstructed feature contrasts, bootstrap stability, and pathway/module enrichment". Since 2026-10-31 the status line reads "not met".
- **The evidence is complete.** No input the review needs is pending:

  | case study | design | primary result | size | power at its stage signal | source |
  |---|---|---|---|---|---|
  | SEA-AD astrocytes `cohort_30` (primary) | n = 80, RNA + ATAC pseudobulk, six covariate columns, rank 5, AUROC 0.636, `df` | `delta` / `angle` / `shape` p 0.501 / 0.449 / 0.331 | **measured** under `df`: pooled 0.029–0.058 per cohort table and statistic ([corrected-permutation report](../../../docs/reports/corrected-permutation-2026-10-30.md)) | **bracketed**: null-rate at n = 80 and AUROC 0.52 / 0.58 ([rank-floor report](../../../docs/reports/rank-floor-2026-10-23.md) §4, every rate 0.010–0.062); no power at n = 200 and AUROC 0.553 / 0.729 ([design-guidance pilot](../../../docs/reports/design-guidance-pilot-2026-10-28.md), δ arm) | `sead_sexdiff` `reports/motco-residual-correction-primary-2026-10-31.md` |
  | SEA-AD astrocytes `cohort_50`; microglia `cohort_30` / `cohort_50` | n = 74; 76 / 63 | no rejection | measured for microglia `cohort_30`'s table; **inferred** for `cohort_50` (κ 1.075 / 1.090) | bracketed as above | same |
  | SEA-AD sensitivity arms (forced ranks, 8 preprocessing variants, latent adjustment, donor exclusion) | 40 analyses | no rejection | inferred; `latent_adjusted` `angle` / `shape` inadmissible (`residualized_basis`) | bracketed | same |
  | ADNI primary | n = 1,501, p180 + NMR, six covariate columns, rank 8, AUROC 0.600, uncorrected (κ ≈ 1.003) | p 0.112 / 0.688 / 0.711 | **measured**: 0.034–0.068 against 0.0695 ([ADNI-design report](../../../docs/reports/adni-design-2026-11-02.md) §2) | **measured**: `non_rejections_uninformative`; `delta` at most 0.530, `angle` / `shape` at the null rate (§3) | `adni_sexdiff` `reports/motco-adni-2026-11-02-2.md` |
  | ADNI sensitivities (latent stage, six export variants, CN / LMCI / AD, forced ranks 2–12) | — | `delta` rejects only at forced rank 6 (p 0.047) and under `nmr_all_columns` (0.045) | inferred | measured only for the primary's construction | same |

- **The attribution path cannot run at the case-study configuration.**
  - `analyze_trajectories` raises "Attribution is not available with covariates: it works from unadjusted cell means" whenever covariates are declared (`src/motco/analysis.py`, the guard before fitting).
  - `analyze_orientation_attribution` accepts covariate-adjusted means only as a precomputed `mean_table`, and then requires `bootstrap_replicates = 0`, so stability is unavailable (`src/motco/stats/attribution.py`).
  - Both case-study primaries declare six covariate columns.
  - So G1's "stable drivers" could not be produced for a claim from either primary, even if one existed.
- **The form.** The [Phase 5 exit review](../../../docs/reports/phase5-exit-review-2026-09-10.md) has the form to follow:
  - the decision first;
  - the inputs, with no new measurement;
  - a criteria table;
  - claimable / withheld / not claimable;
  - each deviation with its mechanism and its kind (method revision or claim revision);
  - the next phase's status;
  - the consequent work in order;
  - what the review does not do.

  It was a docs-only commit with no OpenSpec change. This one has a change because the user proposed it as roadmap item 10.

## Goals / Non-Goals

**Goals:**

- Record one dated decision on Phase 6, with each gate part's status and the evidence for it, so the roadmap stops reading the phase as open.
- Fix what the paper's real-data section may claim, what it must withhold and what it must not claim, each traceable to a committed file.
- Predeclare the criterion that lifts the withheld attribution leg, before the work that would lift it is designed.
- Name the capability gap (attribution with covariates) as a consequent change rather than leave it implicit in the lifting criterion.

**Non-Goals:**

- **Re-reading any case-study result under new rules.** The case studies' predeclared reading rules stand, and the review confirms or cites them.
- **Designing the consequent changes.** Covariate-adjusted attribution and the plasmode each get their own proposal, with thresholds, constructions and compute.
- **Revising the Phase 6 gate's text.** A status is recorded, but no criterion is rewritten after the fact.
- **Ruling on Phase 7.** It stays deferred on the Phase 5 review's §5 grounds, which nothing here touches.

## Decisions

### D1. Phase 6 closes on a size-controlled null, and the attribution leg is withheld (user decision)

The decision is: **CLOSED on a size-controlled null; the attribution leg withheld, lifted by a named demonstration or reported as failed on its own terms.**

- *Why close.* Neither cohort can produce a G1 claim at its design. For ADNI this is measured: no statistic reaches the floor. For SEA-AD it is bracketed: no power at AUROC 0.52–0.73 at n = 80 and n = 200. Keeping the phase open for these cohorts would invite more analyses until one rejects, which the case studies' own predeclared rules forbid. The inference leg is complete: two designs, size measured at both primaries, reproducibility shown across machines (SEA-AD) and by a regenerating `--check` report (ADNI).
- *Why withhold rather than waive.* The roadmap's objective 4 ("identify stable feature and pathway drivers after a global trajectory difference is detected") is the part of the method that has no real-data demonstration. Waiving it would make the paper's real-data section silent on the method's interpretive half. Withholding it keeps the claim open with a named path, as the Phase 5 review did for magnitude specificity: that claim was withheld on 2026-09-10 and lifted by the 2026-10-08 addendum.
- *Alternatives considered* (offered to the user):
  - *close on the null only*: the same closure, with real-data attribution left as a permanent limitation;
  - *hold open until the plasmode has run*: the same work, but Phase 6 would stay open over a leg that needs a MOTCO capability change first;
  - *hold open for a third, powered cohort*: the stage signal limits power, so this means a different biological question, not more AD samples, and a new case study of unknown cost.

### D2. Gate status uses three values, and the vacuous pass is refused

Each gate part is recorded as **met**, **not met** or **unexercised**.
- G2 is *met*.
- G1 is *unexercised*: no claim exists. The review also says that attribution is *unavailable* at the case-study configuration, the D5 gap.
- The gate as a whole is *not met as written*.

The review states explicitly that a literal reading ("each claim" over an empty set) would pass vacuously, and refuses it. The purpose sentence and the conditional Work bullet both presuppose a detected difference, and passing a gate because nothing was found would make the gate meaningless.

- *Alternative:* record G1 as *failed*. That is rejected, because failure would say the attribution was tried on a claim and did not produce stable drivers, and nothing of the kind happened.

### D3. Every evidence statement carries its provenance class

Every size and power statement in the ledger (Context table) is labelled with one of three classes:

- **measured**: a committed decision file at that design and configuration;
- **bracketed**: measured at neighbouring designs whose stage signal or n encloses the cohort's;
- **inferred**: from a mechanism, such as κ's dependence on n and the reduced rank only.

The paper may state a measured or bracketed power statement, naming the class. It may not turn a bracket into a number for the cohort.

- *Why.* SEA-AD's "uninformative" reading rests on a bracket, ADNI's on a measurement, and the two must not be presented as the same kind of evidence.
- *Alternative:* a single "supported / unsupported" label. Rejected because it hides exactly the distinction a reviewer would ask about.

### D4. The cross-study finding becomes a binding reading

The review states the finding the two cohorts share with the design-guidance profiles: **at these designs the held-out stage signal, not n, limits power.**

| run | n | stage AUROC | best power |
|---|---|---|---|
| ADNI | 1,501 | 0.600 | `delta` at most 0.530, at a doubled stage effect |
| design-guidance slice | 80 | 1.000 | `delta` 0.62 at e = 0.1 |
| design-guidance δ arm | 200 | 0.553 / 0.729 | no statistic above 0.09 |

It then adds one binding reading for the paper and any later real-data use of `motco analyze`: *a non-rejection is informative only where a measured operating characteristic at the cohort's own design and stage signal says so; the held-out CV stage AUROC is reported beside the CV-selected rank with every result*. This extends the Phase 6 statement that already requires the rank, without changing it.

- *Alternative:* state the finding descriptively only. Rejected because the Phase 5 review's precedent turns the transferable observable (there, the eigengap) into a reading condition, and the AUROC plays that role here.

### D5. Lifting the attribution leg requires the case-study configuration, so covariate-adjusted attribution comes first

The lifting demonstration must run at the configuration a real claim would be read in: covariates declared, the rank floor, net displacement and, at small n, `df`. Attribution as built cannot do that (Context). The review therefore names two consequent changes in order:

1. **Covariate-adjusted attribution** in `motco analyze` (roadmap item 11). It has three parts:
   - attribution from covariate-adjusted cell means, evaluated at the record's reference row;
   - bootstrap stability that refits the adjustment within each resample on the frozen basis;
   - its own validation that the attribution of an adjusted contrast is stable and unbiased under a confounded covariate.

   It is byte-identical without covariates. Its design is its own proposal.
2. **The plasmode demonstration** (roadmap item 12), which lifts the leg or reports it failed. It is designed in D6.

- *Alternative: a plasmode without covariates.* It is cheaper and runs today. But it would demonstrate attribution of an analysis neither case study runs. Under a confounded covariate the unadjusted contrast is a different estimand: the covariate Type I check's unadjusted arm rejected `shape` at 0.364 on null data. So it could not lift a leg that is withheld for the adjusted analysis. The plasmode's own proposal may still add an unadjusted arm as a descriptive comparison.
- *Alternative: precomputed adjusted means with stability off.* Rejected because G1 asks for *stable* drivers.

### D6. The plasmode runs on SEA-AD astrocyte `cohort_30`, with the construction in MOTCO

What the review fixes now, as the lifting criterion:

- **Data.** The real exported blocks, covariates and design of SEA-AD astrocyte `cohort_30`, the primary cohort:
  - SEA-AD is open data, while ADNI's is under a data-use agreement;
  - its configuration's size under `df` is measured;
  - ADNI's CV rank at the analogue's cap makes its latent geometry the less understood of the two.
- **Construction.** A planted group × stage difference of a declared kind, with a known feature set, added to the real data at declared strengths. The orientation kind comes first, because the Work bullet is conditioned on orientation. Each plasmode unit is analysed by `motco analyze` at the case-study configuration with attribution.
- **Lift criteria.** The leg is lifted when three things hold together:
  - (i) the unplanted construction holds size at the one-sided bound α + 2·√(α(1 − α)/R);
  - (ii) the planted difference is detected at a rate fixed in item 12's proposal;
  - (iii) the planted features are recovered among the attribution's top drivers, with bootstrap sign stability and top-k frequency above thresholds fixed in that proposal before any run.

  Otherwise the leg is reported failed on its own terms, with the failing criterion named.
- **Ownership.** MOTCO owns the planting construction and the decision rule, because it is a method validation like its other checks. It reads `sead_sexdiff`'s export through a documented input contract and commits aggregates only. `sead_sexdiff` supplies the export, and no participant-level data enter MOTCO.
  - *Alternative:* run it entirely in `sead_sexdiff`. That is possible, but it would split MOTCO's validation record across repositories. Item 12's proposal confirms the split.

### D7. The ADNI sensitivity-arm `delta` rejections stay as their rules read them

The review records ADNI's two sensitivity rejections and does not promote them:
- `delta` at forced rank 6 (p 0.047), one of six forced ranks;
- `delta` under `nmr_all_columns` (p 0.045), one of six export variants.

Under `adni_sexdiff`'s predeclared sweep and membership rules they are setting dependence, and no rule counts the arms that reject. The ADNI-design check shows size holds and `delta` is conservative (0.034) at the design, so the rejections are not explained by a liberal null. They are also not a claim, because no multiplicity-controlled rule makes one, and the primary's non-rejection is the result.

- *Alternative:* treat them as a G1 candidate and run attribution on them. Rejected for two reasons. It would be the forking path D1 closes. And the attribution would hit the D5 gap anyway.

### D8. The review sorts the remaining roadmap items

The review classifies every open Next changes item as one of three kinds:

- **owed by Phase 6's withheld leg**: items 11 and 12 (new);
- **owed by a case-study repository**: `adni_sexdiff`'s rule 3 confirmation, citing `adni_design_decision.json`;
- **independent of Phase 6**:
  - item 5, block weighting;
  - item 6, the materiality-rule default;
  - item 7(b), group-specific dispersion;
  - item 8's follow-ups: the `df` default, and the unmeasured `cohort_50` tables, on which no rejection depends;
  - the toy-fixture decision.

The roadmap marks the independent items as not blocking anything Phase 6 owes.

### D9. Report layout

`docs/reports/phase6-exit-review-<date>.md` has these sections:

0. Decision paragraph, then inputs.
1. The exit-gate criteria: G1, G2, the purpose and the six Work bullets, each with its status and evidence.
2. The evidence ledger per case study, from the Context table, with provenance classes.
3. What Phase 6 establishes, withholds, and does not allow to be claimed.
4. Resolution of each deviation, each with its mechanism and kind:
   - the absent claim;
   - the rank-1 degeneracy (method revision, 0.10.0);
   - the residual-df shrinkage (method revision, 0.13.0, re-issue with no status change);
   - the technical axes the covariate-blind basis follows (SEA-AD library quality; ADNI phase) and the inadmissible `stage` adjustment for `angle` / `shape` (claim revision);
   - the ADNI sensitivity rejections (D7);
   - attribution unavailable with covariates (capability gap → item 11);
   - the wider structural analogue at ADNI (claim scope).
5. The cross-study finding and the binding reading (D4).
6. The withheld leg and its lifting criterion (D5, D6).
7. Phase 7: unchanged.
8. Consequent work in order (D8).
9. What the review does not do.

Each number cites the committed file it is read from.

### D10. Roadmap edits are confined to status, items and the deferred list

`docs/roadmap.md` changes in these places:

- the "Not yet established" bullet on real data;
- the Phase 6 exit-gate status line, which gains the decision, the date and the link;
- the Phase 6 Work bullet for attribution, annotated *withheld; item 12*;
- Next changes: item 10 done; items 11 and 12 added; the independence note on items 5, 6, 7(b) and 8;
- "Explicitly deferred", which gains one entry: reading a non-rejection at an unmeasured design or stage signal as evidence that trajectories agree.

Objective 4, the gate's text and every Phase 6 measurement paragraph stay as written. The review is not added to the mkdocs nav, as the Phase 5 review is not. AGENTS.md is not touched: no architecture, default or convention changes.

## Risks / Trade-offs

- **[A number in the review drifts from its source]** → Every number is copied from a committed file (MOTCO decision JSONs and reports; the case-study reports at their committed revisions, recorded in the review's inputs). A task cross-checks each against its source before commit.
- **[The case-study repos change after the review]** (for example, `adni_sexdiff` confirms rule 3) → The review pins each case-study input to a commit hash. Later re-issues do not change the review, which is dated.
- **[Closing reads as declaring the method validated on real data]** → The decision sentence leads with "size-controlled null" and "withheld". The claimable list excludes driver tracing, and the not-claimable list excludes agreement readings.
- **[Item 11 proves harder than expected, and the leg stays withheld indefinitely]** → That is acceptable by design. The paper reports the leg as withheld, with the gap named. The review does not make closing depend on it.
- **[A reviewer of the paper reads the SEA-AD bracket as a measurement]** → D3's provenance classes appear in the ledger and the paper inherits them.
- **[The binding AUROC reading over-generalizes from three regimes]** → It is stated as a reading condition ("informative only where measured"), not as a power model, and the three regimes are listed with their n and stage AUROC.

## Migration Plan

Docs only. The work lands on a branch from `main` (which already merges PR #71) and goes through a PR. Rollback is a revert. No release, no version change, and no fixture regeneration.

## Open Questions

- **The null construction for the plasmode's size criterion (D6 (i)).** Real sex labels may carry a real difference, so "unplanted" data need a null-making step, such as within-stage label permutation that preserves the covariates. Item 12's proposal decides it. It does not change this review, which fixes only that a size criterion exists.
- **Whether item 11's validation needs cluster compute.** Item 11's proposal decides that, after a smoke timing.
