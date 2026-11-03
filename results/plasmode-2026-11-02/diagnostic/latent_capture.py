"""Label-free diagnostic: how much of a planted direction lies inside the rank-5 stage-PLS span.

Fits the case-study latent space at a forced rank of 5 (the primary's CV rank) on the unplanted export, under the
fidelity gate's permuted labels (the latent space does not depend on the labels), and measures, in analysed
coordinates, the share of the norm inside the span of the PLS rotations of (a) the pooled stage displacement and
(b) 200 plant directions drawn as the check draws them (10 features per block, signs, ``sigma_f``). Also reports the
plant's norm at e = 1 against the pooled displacement's.

    python results/plasmode-2026-11-02/diagnostic/latent_capture.py --export-dir <cohort_30> --out <json>
"""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import replace
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "scripts"))
import plasmode_check as check  # noqa: E402

PLANTS = 200
FORCED_RANK = 5


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--export-dir", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    export = check.load_export(args.export_dir)
    labels, _ = check.draw_labels(export.real_labels, export.stages, export.stage_order, check.BASE_SEED, "fidelity", 0)
    metadata = export.metadata.copy()
    metadata[export.group_col] = labels.astype(str)
    settings = {**check.base_settings(export), "attribution": None}
    params = replace(check.analysis_params(settings, None, permutations=0), forced_components=FORCED_RANK)
    _, result = check.analyse(export.blocks, metadata, export, params)
    basis, _ = np.linalg.qr(np.asarray(result.latent_artifacts.model.x_rotations_))
    geometry = export.geometry()
    displacement = geometry["displacement"]
    names = {name: list(export.blocks[name].columns) for name in check.BLOCKS}
    shares, norms = [], []
    for seed in range(PLANTS):
        plant = check.draw_plant(names, seed, check.PLANTED_PER_BLOCK)
        direction = np.zeros_like(displacement)
        offset = 0
        for name in check.BLOCKS:
            indices = np.asarray(plant[name]["indices"], dtype=int)
            direction[offset + indices] = export.sigma()[name][indices] * np.asarray(plant[name]["signs"], dtype=float)
            offset += export.blocks[name].shape[1]
        direction = direction / geometry["scale"]
        norms.append(float(np.linalg.norm(direction)))
        shares.append(float(np.linalg.norm(basis.T @ direction) / np.linalg.norm(direction)))
    payload = {
        "forced_rank": FORCED_RANK,
        "n_features": int(displacement.size),
        "export_files": export.files,
        "pooled_displacement_share_in_span": float(
            np.linalg.norm(basis.T @ displacement) / np.linalg.norm(displacement)
        ),
        "pooled_displacement_norm": float(np.linalg.norm(displacement)),
        "plants": PLANTS,
        "plant_share_in_span": {
            "median": float(np.median(shares)),
            "q05": float(np.quantile(shares, 0.05)),
            "q95": float(np.quantile(shares, 0.95)),
            "random_direction_reference": float(np.sqrt(FORCED_RANK / displacement.size)),
        },
        "plant_norm_at_e1": {"median": float(np.median(norms)), "q05": float(np.quantile(norms, 0.05)),
                             "q95": float(np.quantile(norms, 0.95))},
        "units": "analysed coordinates (pooled per-feature SD, ddof 0); the plant applies to the Male donors only",
    }
    args.out.write_text(json.dumps(payload, indent=1, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=1, sort_keys=True))


if __name__ == "__main__":
    main()
