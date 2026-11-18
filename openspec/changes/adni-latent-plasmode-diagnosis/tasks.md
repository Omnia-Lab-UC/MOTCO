# Tasks

## 1. Driver

- [x] 1.1 Create `scripts/adni_latent_plasmode_diagnosis.py` with `run` and `decide` (design D1, D7): `--records-dir` and `--out-dir` required. `settings.json` binds the parent's five committed files by sha256, re-checked by `run` and `decide`. The driver refuses a changed parent file by name, a record under other settings and a missing or duplicated unit. Verify each refusal in `tests/test_adni_latent_plasmode_diagnosis.py` on hand-built records.
- [x] 1.2 Implement the reproduction (D1): the parent's `calibration_rows` and `capture_rows`, written by its CSV writer, against the committed files byte for byte. Verify a pass on hand-built records whose parent tables were written by the same builders, and `reproduction_failed` after one record is changed.
- [x] 1.3 Implement F1–F4 (D2–D5): the variance components, the delta method, the cluster bootstrap and Fieller's interval. Verify the components against a hand computation, the delta method against its formula, each kind of Fieller interval, and that the component and per-cohort bias columns equal the hand-computed means and the parent's calibration table.
- [x] 1.4 Implement the reading (D6) and the decision file. Verify each verdict from synthetic findings, and that no output holds a per-unit or per-cohort value.

## 2. Run and report

- [ ] 2.1 Run `run` and `decide` on ing from a clean worktree at a pushed revision, on `/home1/tgonza/MOTCO-adni-latent-run/out`, with the output outside the repository. Verify the reproduction, and record the records' sha256 and the revision.
- [ ] 2.2 Commit the outputs under `results/adni-latent-plasmode-diagnosis-2026-11-18/` with `PROVENANCE.txt`. Write `docs/reports/adni-latent-plasmode-diagnosis-2026-11-18.md` with the verdict, the four findings, the stand-ins, what follows and the limitations, and add it to `mkdocs.yml`. Verify that every number the report quotes is in a committed output.
- [ ] 2.3 Add the reading to `docs/measurement-visibility.md` (the latent-plasmode section), update roadmap open item 2, the Done-list entry 23 and the span-limits entry, and add the driver's row to `AGENTS.md`. Verify that `mkdocs build` adds no warnings beyond those on `main`, and that the pre-commit gate passes.

## Workflow follow-up

- Sync the delta specs and archive the change with `adni-latent-plasmode-capture` once both are reviewed.
- Under `information_limited`, the capture waits for a larger cohort or several cohorts; ROSMAP's export is the next candidate.
- Under `resolvable`, a rule version 2 with Fieller's interval is its own change, written before it reads.
