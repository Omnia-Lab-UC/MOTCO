"""Diagnostic behind the 2026-11-02 gate amendment (design D10): why the unadjusted top-20 is not
contaminated at orientation e = 1.0.

For six replicates per (delta, s_z, cell) it runs the check's own unit evaluation at 0 permutations
and B = 0 and reports, on the observed principal contrast: unadjusted - oracle contamination and
projection, oracle truth precision, and the best rank an injected non-driver feature reaches in
the unadjusted |contrast|. Indicative only (six replicates, no permutations, no bootstrap).

Run on 2026-11-02 at commit 9dbac9a (MOTCO 0.14.0), from the repository root, as

    OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
        uv run python results/adjusted-attribution-2026-11-02/calibration/probe/contamination_probe.py \
        2.0,1.0,0.5,0.25 5,10 > contamination_probe.jsonl

``contamination_probe_summary.csv`` holds the per-(delta, s_z, cell) means.
"""

from __future__ import annotations

import json
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[4] / "scripts"))
import adjusted_attribution_check as aa  # noqa: E402


def settings_for(delta: float, s_z: float) -> dict:
    args = aa.build_parser().parse_args(
        ["run", "--out-dir", "unused", "--reps", "6", "--perms", "0", "--bootstrap", "0", "--s-z", str(s_z)]
    )
    settings = aa.build_settings(args)
    for field in ("delta_methyl", "delta_expr", "delta_protein"):
        settings["generator"][field] = float(delta)
    return settings


def probe(job: tuple[float, float, str, float, int]) -> dict:
    from motco.simulations.attribution_diagnostics import derive_truth_driver_features
    from motco.stats.attribution import analyze_orientation_attribution

    delta, s_z, mode, effect, replicate = job
    settings = settings_for(delta, s_z)
    built = aa.build_unit(mode, effect, replicate, settings)
    dataset, result = aa.analyse(built, settings)
    oracle_x, injection = aa.oracle_inputs(built, dataset, result, settings)
    features = result.latent_artifacts.features
    labels = dataset.metadata[["group", "stage"]].astype(str).set_axis(features.index)
    kwargs = dict(groups=list(dataset.groups), stages=list(dataset.stages), bootstrap_replicates=0)
    model = result.latent_artifacts.model
    unadjusted = analyze_orientation_attribution(features, labels, model, **kwargs)
    oracle = analyze_orientation_attribution(oracle_x, labels, model, **kwargs)
    adjusted = result.attribution
    truth = derive_truth_driver_features(built["dataset"], list(features.columns), list(dataset.stages), settings["layers"])
    injected = injection != 0
    contaminant = injected & ~truth.global_mask
    direction = injection / np.linalg.norm(injection)
    out = {"delta": delta, "s_z": s_z, "cell": aa.cell_label(mode, effect), "replicate": replicate,
           "selected_lv": result.record["latent"]["selected_lv"], "cv_auroc": result.record["latent"].get("cv_mean_auroc"),
           "n_truth": int(truth.global_mask.sum()), "n_contaminant": int(contaminant.sum())}
    for name, res in (("unadjusted", unadjusted), ("oracle", oracle), ("adjusted", adjusted)):
        q = res.principal_orientation.observed.directional_contrast
        if q is None:
            out[name] = None
            continue
        order = np.argsort(-np.abs(q), kind="stable")
        top = order[:20]
        rank = np.empty(len(q), dtype=int)
        rank[order] = np.arange(1, len(q) + 1)
        out[name] = {
            "contamination": float(contaminant[top].mean()),
            "truth_precision": float(truth.global_mask[top].sum() / min(20, max(1, int(truth.global_mask.sum())))),
            "projection": float(q @ direction),
            "best_contaminant_rank": int(rank[contaminant].min()) if contaminant.any() else None,
            "q20": float(np.abs(q[order[19]])),
            "max_contaminant": float(np.abs(q[contaminant]).max()) if contaminant.any() else None,
        }
    return out


def main() -> None:
    deltas = [float(value) for value in sys.argv[1].split(",")]
    strengths = [float(value) for value in sys.argv[2].split(",")]
    cells = [("anchor", 0.0), ("orientation", 0.5), ("orientation", 1.0)]
    jobs = [(d, s, m, e, r) for d in deltas for s in strengths for m, e in cells for r in range(6)]
    with ProcessPoolExecutor(max_workers=14) as pool:
        for row in pool.map(probe, jobs):
            print(json.dumps(row), flush=True)


if __name__ == "__main__":
    main()
