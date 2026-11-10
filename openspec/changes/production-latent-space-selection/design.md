# Design

## Context

See proposal.md (Why) for the motivation. The facts below shape the approach.

- **What the comparison settled.** Every label-free or group-aware space either lost in-span power or inflated size; the production stage-PLS span stays. At δ = 2.0 and n = 200 production reads `delta` 0.74 at magnitude e = 0.05; at the weak design (SEA-AD table, δ = 0.11, held-out AUROC 0.617) production rejects nothing above the null at any of its cells (0.01–0.06), and the fitted span keeps 0.10 of the pooled stage program. Its units are reproducible from `latent_space_comparison.py` (`design_generators`, `generator_seed`, `rrpp_seed`, `unit_dataset`, `cell_spec`), and its committed `rate_table.csv` holds per design, cell, arm and statistic the count, rejections and failures at 100 × 199.
- **What the plasmode settled.** On real SEA-AD noise, a sparse random plant outside the stage program is invisible to every trajectory statistic (0.016–0.048 at 3 within-stage SDs), as the span principle predicts. Its export contract (`plasmode_check.load_export`, file hashes), label permutation (`permute_within_stage`, `draw_labels`, seeds from base 2034 by stream, replicate and purpose), analysis-parameter builder and fidelity values (rank 5, AUROC 0.635556, identical for every anchor unit because the space is stage-supervised and label-free) are reusable unchanged. The `paper` stream's replicates 0–99 are a prefix of the committed runs'.
- **The candidate settings** are integration parameters (`block_weighting`, `sparsity`) in the harness and typed fields on `TrajectoryAnalysisParams`, so one arm is one parameter set in either entry point. Under weighting the PLS fits with `scale=False` and `delta` is in the weighted metric (block-weighting design D2); under sparsity the arm's rank and keep fraction come from one CV grid (regularized-stage-pls design D3).
- **Paired units.** The arms share the dataset and, at one RRPP seed with serial draws, the permutation indices, so the information about a difference between arms is in the discordant units; the orthonormal-coordinates check used an exact paired test under a count-plus-ceiling rule for exactly this.
- **Costs.** The comparison ran 4,200 five-arm units in 37.1 core-hours (about 32 s per unit, one PLS CV). The plasmode analysed 1,500 units in 73 core-hours with 999 permutations and B = 200 attribution; a plain CV on 80 × 10,000 is a fraction of that. A sparse CV over the five-step ladder costs about five plain CVs.

## Goals / Non-Goals

**Goals:**
- Attribute any change in size or power to weighting, to sparsity or to their combination, on matched data and draws.
- Read the regime that matters for production (the weak simulated design and real SEA-AD noise) with a plant every arm can see in principle.
- Decide, under a rule fixed before the run, whether a paper grade is worth its cost and for which arm.

**Non-Goals:**
- The paper grade, the default flip, the `mfa` and joint-scope arms, out-of-span questions, a padded simulator (proposal, Non-goals).
- Resolving rates near the bound at R = 100 (SE 0.022 at α): the pilot's bound is 0.0936 and the gain it can credit is about 0.1; the paper grade resolves them.

## Decisions

**D1. The arms are a 2 × 2 factorial, one fit each.**
- `pls_cv`, `weighted`, `sparse_block`, `weighted_sparse_block`. Without the factorial a gain of the combination cannot be attributed, and a loss of one setting could hide behind the other's gain.
- **Why `equal_variance` and the block scope.** `equal_variance` is the sample-independent kind and the multi-block-PLS convention; the block scope is the one that cannot empty a block and is what distinguishes per-block sparsity from a joint budget. Each is the kind its change recommends; the alternatives are follow-ups only on failure.
- **Why every arm refits.** The arms change the fit's input or regularization, so unlike the comparison's `pls_residual_pca` nothing can be built from the production fit. Cost is accepted and measured.
- **Alternative: a design-grid column per arm in the study engine.** `evaluation.integration_params.block_weighting` and `.sparsity` are axes already, so a config could do this at paper grade. The pilot needs the plasmode profile and the paired exact tests, which the engine does not have; the paper grade may use the engine.

**D2. The simulation profile is the comparison's in-span cells plus a gradient.**
- The seven shared cells keep the comparison's names, so `cell_spec`, `generator_seed` and `rrpp_seed` give identical datasets and draws and the production arm reproduces 63 committed rows exactly — the identity check that the harness, the generator and the environment are unchanged. The added cells (`magnitude_e0.25`, `orient_e1.0_phi0`, `shape_e1.0`) are new names with their own RRPP seeds; at the weak design they are where a candidate could first show power, and at δ = 2.0 they add a saturated end to S2 and large plants to S3.
- **Why no φ > 0 or sparse-support cells.** Every arm is stage-supervised; the comparison showed those plants keep ≤ 0.075 of their energy in any such span. They would measure nothing and cost a third of the units.
- **Why 100 × 199.** It is the parent's grade and makes the reproduction exact. The pilot decides whether a paper grade is run, not what production is.

**D3. The plasmode profile plants in the stage span on real noise.**
- **The plant.** In the analysis's standardized units (pooled mean and SD over all donors, the sample scaler `analyze_trajectories` fits), let m̄_s be the pooled stage-s mean over both groups and m̄ the mean of the m̄_s. Every row labelled `Male` at stage s gains e · (m̄_s − m̄) per feature, mapped back through the block's pooled SD. The planted group's configuration is then the pooled one stretched by 1 + e about the centroid: the stage-span magnitude edit (`edit_stage_configuration`'s rule) with the sample configuration for the population one. It is in the fitted span to the extent the span captures the pooled stage program, which is what the arms differ on.
- **Why magnitude only.** It is the one edit definable from the pooled sample configuration without a second principal direction (orientation) or an interior stage's off-chord direction (shape), and `delta` is the statistic with power at the plasmode's n. The simulation profile carries the other two.
- **Why the labels are permuted first.** As in the plasmode: the real labels are never analysed, so the anchor is a null and no record can hold a donor id or label vector. The plant goes on the permuted `Male` set, which makes B − A a stage-dependent shift and A the untouched pooled configuration up to sampling.
- **Fidelity as the identity check.** The anchor's production arm has one label-free latent space whatever the permutation; its rank and AUROC equal the plasmode's committed fidelity values or the environment has drifted. This replaces a 500-replicate reproduction of `plasmode_summary.csv`, which R = 100 cannot provide and which would cost 400 further units per arm.
- **Why 999 permutations here.** RRPP is cheap at n = 80 and rank ≈ 5, and the plasmode's committed rates are at 999; the cost is in the fits.

**D4. The decision pairs where it can and bounds where it must.**
- **S1 and S3 are bounds, not paired tests.** Size and specificity are properties of an arm, not of its difference from production: a candidate that is liberal where production is conservative must fail on its own rate. The bound is the comparison's (α + 2·SE at R = 100), the rule the parent's cells were read under.
- **S2 and G are exact paired tests.** The arms share data and draws, so P(Binomial(b + c, ½) ≥ b) on the discordant counts is the right test and is defined when b + c is small (p = 1 at b = c = 0). S2 uses the orthonormal check's count-plus-ceiling rule over 16 cells; a per-cell 2-SE rule would false-fail with probability about 1 − 0.977¹⁶ ≈ 0.31.
- **G asks for one Bonferroni-protected gain.** At R = 100 the smallest credible paired gain is large (about 10 discordant units in one direction), and the 11 gain cells are at the regime where production reads the null rate, so one protected gain is strong evidence and a pooled 2-SE rule would credit a diffuse drift. The pooled gain is reported for the reader; it does not decide.
- **Order and preference.** Size before power before specificity before gain: a candidate liberal at the null is not a candidate. Among qualifying arms the simplest wins (`weighted`: one scalar per block; `sparse_block`: a CV grid; the combination last), because the paper grade and the default flip cost the same whichever is named and a simpler production rule is easier to carry.
- **Alternative: a gain criterion on capture.** `pooled_stage` capture at the weak design is the mechanism the arms are meant to move, but it is not available on real data and a capture gain without a rate gain buys production nothing. It is descriptive.

**D5. Execution follows the recent validation drivers.**
- One script with `run --profile` / `merge --profile` / `decide`; `settings.json` bound by hash with the parents' sha256 (comparison `settings.json`, `calibration/calibration.json`, `rate_table.csv`; plasmode `settings.json` and the export's five files); round-robin shards; resume by complete unit key; records as JSONL; an sbatch array with `PLS_PROFILE`, BLAS pinned to one thread, run on `ing` from a clean worktree at a pushed revision that carries 0.22.0.
- What gets committed: aggregates, the decision, provenance (`PROVENANCE.txt`, `sacct.txt`, package versions) and the report under `results/production-latent-space-<date>/`; records stay on the cluster.

## Risks / Trade-offs

- **[The simulation's 498 features are not the p ≫ n regime where sparsity pays]** → By design the plasmode carries that regime (10,000 features, n = 80). A `no_gain` in simulation with a gain on the plasmode is a readable outcome: G is one criterion over both profiles.
- **[A sparse arm selects an unstable fraction on the weak design and its rank jumps]** → The cell table reports the rank and fraction distributions; instability is visible and is itself a reason not to propose the arm.
- **[Weighted `delta` is in another metric]** → Only p-values and rejections are compared; the metric note is in the decision and the report.
- **[Compute overrun]** → A local smoke (2 replicates per cell, both profiles) measures per-arm time before launch. Levers: a three-step ladder for the sparse arms, 50 plasmode replicates for the plant cells. The projected core-hours go to the user with a cap.
- **[The export changes]** → The hashes are the plasmode's committed ones; a differing file is refused by name.
- **[Environment drift since 2026-11-04]** → The simulation reproduction fails and gates the reading; run on the same `ing` environment and record versions.
- **[The real-noise plant partly follows the sex-permuted labels' sampling noise]** → The plant is defined on the pooled stage means over both groups, so the permutation of sex does not enter it; the anchor's size reading covers the rest.

## Migration Plan

Nothing migrates; no package output changes.
- **On `propose_paper_grade`:** a separate change runs the named arm at 500 × 999 (the stage-span paper grade's cells; the plasmode at its grade), and on its verdict another flips the harness and `motco analyze` defaults with pins and the SEA-AD re-issue.
- **On `keep_production`:** both settings stay opt-in and documented as measured; the measurement-visibility page records what they changed and did not.

## Open Questions

- The core-hour cap and whether the ladder or the plasmode replicate count is shortened are set from the smoke timing, at the approval gate.
