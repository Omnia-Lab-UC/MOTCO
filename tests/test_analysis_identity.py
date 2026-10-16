"""No-covariate byte identity of the real-data analysis.

``tests/data/analysis_no_covariate_identity.json`` was captured with the 0.7.0
``analyze_trajectories`` (before covariate support existed) by running this
module as a script::

    uv run python -m tests.test_analysis_identity

It holds the record without ``measurement.runtime_seconds`` and ``software``
(both vary by run or environment) and SHA-256 digests of the three CSV outputs.
The test requires an analysis without covariates to reproduce it exactly.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from motco.analysis import (
    OmicsBlock,
    TrajectoryAnalysisParams,
    TrajectoryDataset,
    analyze_trajectories,
    write_analysis_outputs,
)

FIXTURE = Path(__file__).parent / "data" / "analysis_no_covariate_identity.json"
STAGES = ["Low", "Intermediate", "High"]
CSV_OUTPUTS = ("latent_scores.csv", "pls_loadings.csv", "stage_means.csv")


def _identity_dataset() -> TrajectoryDataset:
    rng = np.random.default_rng(20261016)
    rows, rna, atac = [], [], []
    for stage_index, stage in enumerate(STAGES):
        for group_index, group in enumerate(("F", "M")):
            for replicate in range(7):
                rows.append((f"{group}_{stage}_{replicate}", group, stage))
                shift = float(stage_index) * (1.0 + 0.4 * group_index)
                rna.append(shift + rng.normal(scale=0.7, size=12) + np.linspace(0, group_index * stage_index, 12))
                atac.append(0.5 * shift + rng.normal(scale=0.7, size=8))
    ids = [row[0] for row in rows]
    metadata = pd.DataFrame({"donor": ids, "Sex": [r[1] for r in rows], "ADNC": [r[2] for r in rows]})
    return TrajectoryDataset(
        blocks=[
            OmicsBlock("rna", pd.DataFrame(np.asarray(rna), index=ids, columns=[f"g{j}" for j in range(12)])),
            OmicsBlock("atac", pd.DataFrame(np.asarray(atac), index=ids, columns=[f"p{j}" for j in range(8)])),
        ],
        metadata=metadata,
        sample_col="donor",
        group_col="Sex",
        stage_col="ADNC",
        stage_order=STAGES,
    )


def _identity_params() -> TrajectoryAnalysisParams:
    return TrajectoryAnalysisParams(
        permutations=99, seed=7, cv2_splits=3, cv1_splits=2, n_repeats=2, max_components=5
    )


def capture_identity(out_dir: Path) -> dict[str, Any]:
    """Run the fixed analysis, write its artifact to ``out_dir`` and return the identity payload."""

    result = analyze_trajectories(_identity_dataset(), _identity_params())
    paths = write_analysis_outputs(result, out_dir, overwrite=True)
    record = json.loads(paths["result"].read_text(encoding="utf-8"))
    record.pop("software")
    record["measurement"].pop("runtime_seconds")
    digests = {
        name: hashlib.sha256((Path(out_dir) / name).read_bytes()).hexdigest() for name in CSV_OUTPUTS
    }
    return {"record": record, "csv_sha256": digests}


def test_no_covariate_analysis_is_byte_identical_to_0_7_0(tmp_path: Path) -> None:
    expected = json.loads(FIXTURE.read_text(encoding="utf-8"))
    observed = capture_identity(tmp_path)
    assert "covariates" not in observed["record"]
    assert observed["record"] == expected["record"]
    assert observed["csv_sha256"] == expected["csv_sha256"]


if __name__ == "__main__":
    import tempfile

    with tempfile.TemporaryDirectory() as scratch:
        payload = capture_identity(Path(scratch))
    FIXTURE.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"wrote {FIXTURE}")
