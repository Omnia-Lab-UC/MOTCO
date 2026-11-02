# phase6-exit-review

## Why

Phase 6's exit gate ("collaborators can trace each global trajectory claim to stable feature/module drivers and reproduce the result from documented inputs and configuration") has been recorded as *not met* since 2026-10-31, because neither case study has a claim to trace:

- **SEA-AD** (`sead_sexdiff`, n = 80 / 74 astrocytes and 76 / 63 microglia, MOTCO 0.13.0 with `--residual-correction df`) rejects no statistic in any of its 44 analyses.
- **ADNI** (`adni_sexdiff`, n = 1,501, MOTCO 0.10.0) rejects nothing in its primary. `delta` rejects only at forced rank 6 (p 0.047) and under the `nmr_all_columns` export variant (0.045). Its predeclared rules read both as setting dependence.

`motco analyze --attribution` has therefore never run on real data. Roadmap Next changes item 10 asks for a dated exit review, in the form of the [Phase 5 exit review](../../../docs/reports/phase5-exit-review-2026-09-10.md), that decides how Phase 6 closes. The last measurement the review waited for is item 9's ADNI-design check: size holds, and no statistic reaches the 0.80 floor (`non_rejections_uninformative`, 2026-11-02, PR #71, merged). Without a dated decision three things go wrong:
- Phase 6 reads as open indefinitely.
- The paper's real-data section has no recorded basis for what it may claim.
- Any further analysis of the two cohorts aimed at producing a claim would be a forking path.

## What Changes

- **A dated exit review**, `docs/reports/phase6-exit-review-<date>.md`, in the Phase 5 review's form. The decision comes first, then the inputs, then a criteria table, then what is claimable, withheld and not claimable, the resolution of each deviation, the consequent work in order, and what the review does not do. It makes no new measurement. Every number it states is read from a committed MOTCO report or decision file, or from the case-study repositories' committed reports, and cited by path.
- **The decision**, taken with the user before this proposal: **Phase 6 is CLOSED on a size-controlled null, with the attribution leg withheld.**
  - The gate is recorded as *not met as written*. Its reproducibility criterion is met. Its driver-tracing criterion is *unexercised*, because there is no claim to trace; that is different from failed.
  - The phase's purpose is demonstrated up to inference: an end-to-end, reproducible, size-controlled analysis of two real cohorts whose non-rejections are measured to be uninformative.
  - The attribution leg is withheld, not waived. A named demonstration lifts it, or reports it as failed on its own terms. This mirrors the Phase 5 review's withheld magnitude-specificity claim and its 2026-10-08 addendum.
- **An evidence ledger per case study.** Each design, primary result and sensitivity result is recorded with:
  - where its size comes from: measured (SEA-AD `cohort_30` under `df`; ADNI) or inferred (SEA-AD `cohort_50`, the forced ranks, the variants);
  - where its power statement comes from: measured (ADNI) or bracketed (SEA-AD, by the rank-floor check at n = 80 and the design-guidance δ arm at n = 200);
  - its binding reading.
- **The cross-study finding stated for the paper.** At these designs the held-out stage signal limits power more than n does. ADNI at n = 1,501 and AUROC 0.60 has no useful power. The design-guidance slice at n = 80 and AUROC 1.000 reaches `delta` 0.62 at e = 0.1.
- **A binding reading for any later real-data use.** A non-rejection is informative only where a measured operating characteristic at the cohort's own design and stage signal says so. The held-out CV stage AUROC is reported beside the CV-selected rank with every result.
- **The lifting criterion, predeclared, and the gap it exposes.** `analyze_trajectories` refuses attribution whenever covariates are declared (`src/motco/analysis.py`, "Attribution is not available with covariates"). `analyze_orientation_attribution` accepts precomputed adjusted means only with bootstrap stability disabled. Both case-study primaries are covariate-adjusted, so no claim from either could be traced to *stable* drivers as the analysis is built today. The review names two consequent changes, each with its own proposal:
  - **covariate-adjusted attribution** with bootstrap stability in `motco analyze` (MOTCO);
  - **a plasmode demonstration** on SEA-AD's real astrocyte blocks: a planted group × stage difference with known features is analysed at the case-study configuration with attribution. It lifts the withheld leg if the planted difference is detected and the planted features are recovered with bootstrap stability while the unplanted data hold size.
- **Roadmap.** `docs/roadmap.md` changes in six places:
  - the Phase 6 exit-gate status line gains the decision and a link;
  - the Phase 6 Work bullet for attribution is annotated;
  - item 10 is marked done;
  - the two consequent changes are added as items 11 and 12;
  - the remaining items (5, 6, 7(b), 8's follow-ups) are marked as not owed by Phase 6;
  - "Explicitly deferred" gains "Reading a non-rejection at an unmeasured design or stage signal as evidence that trajectories agree".

  The Phase 6 gate's text and every predeclared artefact stay as written.
- **Out of scope.** The review makes no change under `src/motco/` or `tests/`, runs no measurement and creates no release. It does not edit `sead_sexdiff` or `adni_sexdiff`. ADNI's reading rule 3 confirmation stays that repository's own change, and the review lists it as consequent work. Phase 7 stays deferred on the Phase 5 review's grounds.

## Capabilities

### New Capabilities

None.

### Modified Capabilities

None. The review is documentation: it records a decision and changes no system behaviour, so the change sets `skip_specs: true`, as `phase5-paper-grade-run` did. The two consequent changes will each declare their own capabilities: covariate-adjusted attribution would modify `real-data-trajectory-analysis` and `orientation-driver-attribution`, and the plasmode would be a new validation capability.

## Impact

- **Docs:** new `docs/reports/phase6-exit-review-<date>.md`; `docs/roadmap.md` (the Phase 6 section's status line and Work bullet, Next changes items 10–12, the remaining items' status, "Explicitly deferred").
- **Code and tests:** none. The pre-commit gate still runs before the commit.
- **Downstream:**
  - `sead_sexdiff` and `adni_sexdiff` cite the review in their reports' limitations when they next re-issue; nothing in them changes now.
  - `adni_sexdiff`'s rule 3 confirmation (citing `adni_design_decision.json`) stays pending there.
  - The paper's real-data section is bounded by the review's claimable and withheld lists.
- **Compute:** none.
