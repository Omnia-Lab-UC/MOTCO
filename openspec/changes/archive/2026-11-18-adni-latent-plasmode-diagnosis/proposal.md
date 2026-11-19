# Proposal

## Why

The ADNI latent plasmode (`adni-latent-plasmode-capture`, 2026-11-18) read `gate_failed`. In its calibration, the bias of `pls_cv`'s total capture exceeded ±0.05 at five of eight strength × point cells. Recorded with the gate: the program was not present at 80 participants, and no point resolved (half-widths 1.13 to 3.59 against 0.15). Its rule names the consequence: a diagnosis of the estimator comes before any reading.

The run's records already hold what a diagnosis needs. For every unit they hold the numerator and denominator of every reading and, in the calibration, their truths. So the questions can be answered without a new unit:

- Is the bias in the estimator's parts, or in the ratio?
- Does the spread shrink with more units, or is it a property of the cohort?
- Could any rule bound ADNI's capture at all?

An exploratory pass on the same records (`notes.md`; not a result) found the parts unbiased, the spread almost all between cohorts, and ADNI's program energy too close to zero for a bounded interval. This change records that diagnosis under a driver, with its bindings and reproduction, so that the roadmap can state the verdict.

## What Changes

- **A driver, and no package change.** `scripts/adni_latent_plasmode_diagnosis.py` (`run` | `decide`).
    - It reads the latent-plasmode run's merged records where they are (on ing), and binds the parent's committed files by sha256.
    - It refuses records written under other settings, or with a missing unit.
    - It reproduces the committed calibration and capture tables byte for byte from the records before any reading.
- **Four findings on the designated reading** (`pls_cv`, total capture), recorded for every space and stage contrast:
    - F1: the parts are unbiased.
    - F2: the ratio of sums over all cohorts is on its truth.
    - F3: where the spread lives (between and within cohorts).
    - F4: whether ADNI's own 95 % Fieller interval for the capture is bounded, under four stand-ins for its between-cohort covariance.
- **A reading rule** (version 1): `reproduction_failed`, `estimator_biased`, `resolvable`, `units_limited` or `information_limited`, with what each means for a rule version 2. It is written after the exploratory pass and says so; it does not change the parent's verdict.
- **Report and records.** A dated report, aggregates under `results/`, the reading on the measurement-visibility page and the verdict in the roadmap.

## Capabilities

### New Capabilities

- `stage-program-capture-diagnosis`: the diagnosis of a capture estimator's calibration failure from its own records, with its parts, its aggregate, its spread's location and the cohort's own interval.

### Modified Capabilities

- `latent-space-metric-compatibility`: the measurement-visibility page documents the diagnosis and what it means for the real-noise capture.

## Impact

- **Code:** `scripts/adni_latent_plasmode_diagnosis.py` and `tests/test_adni_latent_plasmode_diagnosis.py`. It imports `adni_latent_plasmode_check` unchanged.
- **Package, API, CLI, dependencies:** none. No version bump and no changelog entry.
- **Data:** the run's records stay on ing, and their sha256 is recorded. No export is read. Outputs hold aggregates only.
- **Compute:** about a minute on one core.
- **Docs:** `docs/reports/adni-latent-plasmode-diagnosis-2026-11-18.md`, `docs/measurement-visibility.md`, `docs/roadmap.md`, `mkdocs.yml`, and one row in `AGENTS.md` (Validation drivers).
