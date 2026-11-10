"""Pack the mixOmics sPLS-DA reference fixture for ``SparsePLSDA``.

Run from the repository root after ``tests/data/reference/splsda_fixture.R``::

    Rscript tests/data/reference/splsda_fixture.R build/splsda_fixture
    uv run python -m tests.splsda_fixture_capture build/splsda_fixture

Writes ``tests/data/splsda_reference_fixture.npz``: the inputs (``X``, the
one-hot ``Y`` mixOmics fitted, the class labels), mixOmics's X weights
(``loadings$X``) and X scores (``variates$X``), the rank and the per-component
budget, and the provenance (R and mixOmics versions, the R script's sha256,
mixOmics's iteration counts). ``tests/test_sparse_pls.py`` reads it without R.
"""

from __future__ import annotations

import hashlib
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "data" / "reference" / "splsda_fixture.R"
FIXTURE = ROOT / "data" / "splsda_reference_fixture.npz"


def _provenance(path: Path) -> dict[str, str]:
    entries = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        key, _, value = line.partition("=")
        entries[key.strip()] = value.strip()
    return entries


def capture(export_dir: Path) -> Path:
    provenance = _provenance(export_dir / "provenance.txt")
    np.savez_compressed(
        FIXTURE,
        X=pd.read_csv(export_dir / "X.csv").to_numpy(dtype=float),
        Y=pd.read_csv(export_dir / "Y.csv").to_numpy(dtype=float),
        classes=pd.read_csv(export_dir / "classes.csv")["class"].to_numpy(dtype=str),
        weights=pd.read_csv(export_dir / "weights.csv").to_numpy(dtype=float),
        scores=pd.read_csv(export_dir / "scores.csv").to_numpy(dtype=float),
        keep=np.int64(provenance["keep"]),
        ncomp=np.int64(provenance["ncomp"]),
        r_version=np.str_(provenance["r_version"]),
        mixomics_version=np.str_(provenance["mixomics_version"]),
        mixomics_iterations=np.array([int(v) for v in provenance["iterations"].split(",")], dtype=np.int64),
        script_sha256=np.str_(hashlib.sha256(SCRIPT.read_bytes()).hexdigest()),
    )
    return FIXTURE


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("usage: python -m tests.splsda_fixture_capture <export_dir>")
    print(f"wrote {capture(Path(sys.argv[1]))}")
