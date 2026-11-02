# Tasks

## 1. Preconditions and pinned inputs

- [x] 1.1 Branch `docs/phase6-exit-review` from `main` at or after `0d72692` (PR #71 merged), leaving the untracked `.idea/` and `results/rank-floor-2026-10-23/*` files alone; verify `git merge-base --is-ancestor 0d72692 HEAD` succeeds and `results/adni-design-2026-11-02/run/adni_design_decision.json` reads `"verdict": "non_rejections_uninformative"` with `"indicative": false`.
- [x] 1.2 Pin the case-study inputs (design D9, Risks): record `git -C ~/Documents/Repos/sead_sexdiff rev-parse HEAD` and `git -C ~/Documents/Repos/adni_sexdiff rev-parse HEAD` with each tree's tracked-file status, and the report paths the review reads (`reports/motco-residual-correction-primary-2026-10-31.md`; `reports/motco-adni-2026-11-02-2.md`); verify both reports exist at those revisions (`git -C <repo> cat-file -e <rev>:<path>`).
- [x] 1.3 Confirm the attribution gap the review names (design D5) on the branch's tree: `uv run pytest tests/test_analysis.py::test_attribution_with_covariates_is_rejected_before_fitting -q` passes, and the `bootstrap_replicates must be zero when mean_table contains precomputed means` guard is present in `src/motco/stats/attribution.py`; record both locations (file:line) for the review's §4.

## 2. Write the exit review

- [x] 2.1 Create `docs/reports/phase6-exit-review-<date>.md` (the date it is written) with the decision paragraph verbatim from design D1 and an Inputs paragraph listing every MOTCO report and decision file read and the two case-study reports at the hashes from 1.2, stating that no new measurement was made; verify the decision sentence contains "size-controlled null", "not met as written" and "withheld".
- [x] 2.2 Write §1 (the exit-gate criteria): G1 `unexercised` (and attribution `unavailable` at the case-study configuration), G2 `met`, the purpose `demonstrated up to inference`, and one row per Phase 6 Work bullet, each with its evidence path; include the refusal of the vacuous pass (design D2); verify every row cites a committed file.
- [x] 2.3 Write §2 (the evidence ledger) from design's Context table, one row per case-study analysis group, every size and power entry labelled `measured`, `bracketed` or `inferred` (design D3); verify each label against its source: the corrected-permutation report and the ADNI-design decision for the measured statements, rank-floor §4 and the design-guidance pilot's δ arm for the bracket, κ for `cohort_50`.
- [x] 2.4 Write §3 (claimable / withheld / not claimable). Claimable: end-to-end reproducible analysis of two designs through one entry point; size at both primaries; the null result stated with its uninformativeness; the design finding (D4); the configuration lessons. Withheld: driver tracing on real data. Not claimable: any agreement reading, sensitivity-arm rejections as results (D7), a numeric power for SEA-AD, and metabolite- or feature-level claims from the structural analogue. Verify no item in "claimable" depends on an `inferred` ledger entry.
- [x] 2.5 Write §4 (resolution of each deviation), the seven deviations of design D9 §4, each with its mechanism, its kind (method revision / claim revision / capability gap / claim scope) and its source report; verify the two ADNI sensitivity p-values (0.047, 0.045) and the SEA-AD positive-control outcome (E3 fails on the primary basis in `cohort_30`, holds on `latent_adjusted`) match the case-study reports at the pinned hashes.
- [x] 2.6 Write §5 (the cross-study finding and the binding reading, design D4) with the three regimes' n, stage AUROC and best power, and the reading sentence verbatim from D4; verify each regime's numbers against the ADNI-design report §3, the design-guidance slice report and the design-guidance pilot's δ arm.
- [x] 2.7 Write §6 (the withheld leg and its lifting criterion): the gap (1.3's locations), item 11 (covariate-adjusted attribution with resampled-adjustment bootstrap, byte-identical without covariates) and item 12 (the plasmode on SEA-AD astrocyte `cohort_30`, the three lift criteria (i)–(iii) of design D6, the ownership split, and "reported failed on its own terms"); verify the thresholds are explicitly left to item 12's proposal and fixed there before any run.
- [x] 2.8 Write §7 (Phase 7 unchanged), §8 (consequent work in order, design D8: items 11 → 12; `adni_sexdiff` rule 3; then the independent items 5, 6, 7(b), 8's follow-ups and the toy fixture, marked non-blocking) and §9 (what the review does not do: no measurement, no code, no edit to any predeclared artefact or to the gate's text, no re-reading of case-study results, no further analysis of the two cohorts read as Phase 6 evidence); verify every open Next changes item appears exactly once in §8.
- [x] 2.9 Cross-check every number in the review against its cited source in one pass, using a scratch table of number → file → location kept outside the repo; verify that no number lacks a citation and that every value matches its source to the printed digits.

## 3. Roadmap

- [x] 3.1 Update `docs/roadmap.md`: the Phase 6 exit-gate status line (decision, date, link to the review), and the "Not yet established" bullet on real data (the closure, the withheld leg, items 11–12); verify `grep -n "not met" docs/roadmap.md` shows no Phase 6 status line still reading open.
- [x] 3.2 Annotate the Phase 6 Work bullet "For significant orientation, …" with *withheld — no claim; unavailable with covariates; items 11–12*, leaving the bullet's text and the gate's text unchanged; verify with `git diff` that only the annotation is added to those lines.
- [x] 3.3 In Next changes, mark item 10 done (date, link, decision in one sentence). Add item 11 (covariate-adjusted attribution) and item 12 (plasmode demonstration), each stating its scope boundary from design D5 / D6 and "its own proposal". Add the independence note to items 5, 6, 7(b) and 8's follow-ups. Verify the items run 1–12 without gaps and that each new item links the review.
- [x] 3.4 Add the "Explicitly deferred" entry "Reading a non-rejection at an unmeasured design or stage signal as evidence that trajectories agree"; verify it is the only addition to that list.
- [x] 3.5 Verify every relative link in the new review and in the edited roadmap lines resolves to an existing file (a short link-check over `](...)` targets), and that `uv run mkdocs build` reports no warning that names the review or `roadmap.md` beyond those already present on `main`.

## 4. Integration checks and commit

- [x] 4.1 Run `openspec validate phase6-exit-review --strict` and verify it passes with specs skipped.
- [x] 4.2 Run the pre-commit gate (`uv run ruff check src/ tests/ && uv run mypy src/motco/ && MOTCO_TEST_PERMS=99 uv run pytest tests/ -m "not slow" --tb=short`) and verify all three pass; nothing under `src/` or `tests/` has changed (`git diff --stat main -- src tests` is empty).
- [x] 4.3 Commit as `docs: Phase 6 exit review — closed on a size-controlled null, attribution leg withheld` with the change's artifacts, push the branch and open a PR (problem statement, the decision, the files changed, the commands of 4.1–4.2, the OpenSpec change); verify `git status` shows no stray tracked changes and the PR lists only `docs/` and `openspec/changes/phase6-exit-review/`.

## Workflow follow-up

- After the PR is merged (the merge is the user's), archive the change (`openspec archive phase6-exit-review`; specs are skipped) and verify the archived directory.
- `adni_sexdiff`: confirm reading rule 3, citing `adni_design_decision.json`, in that repository's own change.
- Propose item 11 (covariate-adjusted attribution) and then item 12 (the plasmode), each as its own OpenSpec change.
