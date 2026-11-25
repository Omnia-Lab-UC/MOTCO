from __future__ import annotations

import os
import sys
from pathlib import Path

# One BLAS thread per test process, set before numpy loads: the suite's matrices are
# small, a thread per core only spins, and pytest-xdist runs one process per core.
# An explicit setting in the environment wins.
for _variable in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_variable, "1")

import pandas as pd  # noqa: E402
import pytest  # noqa: E402

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = PROJECT_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))


@pytest.fixture(scope="session")
def data_dir() -> Path:
    return PROJECT_ROOT / "tests" / "data"


@pytest.fixture(scope="session")
def group_candidates() -> tuple[str, ...]:
    # Provided simulations use either "taxa" or "tax" as the group column
    return ("taxa", "tax")


@pytest.fixture(scope="session")
def level_col() -> str:
    # In provided simulation: "Inv" denotes the levels
    return "Inv"


@pytest.fixture(params=["evo_649_sm_example1.csv", "evo_649_sm_example2.csv"], scope="session")
def example_df(request: pytest.FixtureRequest, data_dir: Path) -> pd.DataFrame:
    """Load one of the example CSV datasets as a DataFrame.

    Columns:
    - group column: "taxa"
    - level column: "Inv"
    - remaining numeric columns are features
    """
    csv_path = data_dir / str(request.param)
    df = pd.read_csv(csv_path)
    return df


@pytest.fixture()
def group_col(example_df: pd.DataFrame, group_candidates: tuple[str, ...]) -> str:
    for cand in group_candidates:
        if cand in example_df.columns:
            return cand
    raise AssertionError(
        f"No group column found. Tried {group_candidates}, got {list(example_df.columns)}"
    )
