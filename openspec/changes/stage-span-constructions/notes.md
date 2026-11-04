# Run notes — stage-span-constructions

## Implementation decisions taken while applying

- **Span-capture shares are energy shares.** The delta spec defines `group_difference` as the share of
  Σ_s‖d_s‖² inside the span and `random_share` as rank / n_features, so both are squared-norm (energy) shares.
  The plasmode's `latent_capture.py` reported a norm share ‖Qᵀd‖/‖d‖ (2.16 %, random reference √(5/10,000)).
  For one direction the two are related by a square: `tests/test_span_capture.py` checks
  `span_share(d, Q) == (plasmode share)²` on a fixed matrix, and `docs/measurement-visibility.md` quotes the
  plasmode figure in both units (2.2 % of the norm, 0.05 % of the energy). R2's 0.90 floor is on the energy share.
- **Byte identity is pinned on the BLAS-free parts.** `numpy`'s `multivariate_normal` takes an SVD, so the last bits
  of the sampled matrices depend on the BLAS thread count (the digest differed between one thread and the default).
  `tests/test_stage_span_construction.py` pins, against 85b74bb, the digest of the indicators, metadata and
  population means plus the cell ids and signatures, and compares the sampled matrices in-process. The full sampled
  matrices were checked byte-identical to 85b74bb at one thread and at the default threading.
- **Committed reports regenerate byte-identically.** The six runs with local merged records
  (`phase6-small-n-pilot`, `phase6-small-n-magnitude`, `phase6-small-n`, `phase5-2026-09-10`,
  `net-displacement-type-i/{phase5,phase6}`) were re-reported with 85b74bb and with this change: every file
  (PNGs included) is identical between the two, all 83 committed CSV/JSON files match, and no `span_capture.csv`
  is written.
- **Unit count.** The four profiles enumerate 92 cells = 9,200 units at 100 replicates (14 + 32 per sibling pair
  member). The proposal's 8,400 left out the 8 Type I baseline cells (2 per profile; the null Type I cell is the
  same cell in all four profiles).
- **The prototype was not on disk.** `scripts/latent_retrofit_check.py` was written from the spec: it reuses
  `edit_stage_configuration` on the latent stage means, so the latent edit and the population plant are the same
  geometric operation; edit sizes default to magnitude 0.25, orientation 1/3 (30°), shape 0.5, translation 0.5.

## 10.1 Local smoke (2026-11-04)

Four scratch copies of the profiles at `n_replicates` 2 and 199 → 19 permutations, 6 shards each, 7 worker
processes on the workstation (8 cores, BLAS pinned to one thread), then `merge` and `report` per profile.

- **Units:** 184 / 184 completed (28 + 64 + 28 + 64), 0 failures. Every record carries `realized_geometry.span_capture`
  and the `population_latent` checkpoint. The anchors are `construction`-free (`surgery`/`none`) in all four profiles;
  the stage-span profiles' effect cells and Type I translation record `construction: "stage_span"`.
- **Reports:** all four write `span_capture.csv` beside the usual tables; `anchor_decision` reads `size_holds` at the
  three columns (2 replicates, not informative).
- **Reading:** `scripts/stage_span_pilot_reading.py` runs end to end on the smoke reports. At 2 replicates it reads
  `revise_construction` on R2 alone: stage-span orientation at e = 1.0 has median `group_difference` 0.83 at
  n = 80 and 0.84 at the SEA-AD table. Not a finding at R = 2, but the cell to watch: a 90° rotation in the
  population plane moves the plant furthest from the sample's fitted span at small n.
- **Timing:** median 15.7 s per unit (mean 16.8, max 20.9) at 19 permutations under 7-way SMT contention on the
  laptop; 3,090 core-seconds for the 184 units.

**Projection.** The design-guidance pilot ran the same generator and measurement at n = 80 / 200 at 199
permutations: 10.7 / 13.2 s per unit in its ing smoke and 26 core-seconds per unit over its run (contention
included). Span capture adds one `transform` and one QR per unit and the stage-span plant one population-means
pass and one `estimate_difference` on 6 rows — negligible against the double CV. At 26 core-seconds per unit the
9,200 units need about **66 core-hours**; with a 20 % margin, **80 core-hours**. Shards sized at 30 s per unit:
8 shards per magnitude profile (175 units, about 1.5 h) and 16 per geometry profile (200 units, about 1.7 h),
48 array tasks in all.

## 10.2 Cluster run

(filled when the run completes)
