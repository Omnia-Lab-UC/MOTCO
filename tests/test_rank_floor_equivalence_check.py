"""``scripts/rank_floor_equivalence_check.py``: the summary, the JSONL conversion and the outputs."""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import pytest


def _script(name: str):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
    try:
        return __import__(name)
    finally:
        sys.path.pop(0)


check = _script("rank_floor_equivalence_check")


def _row(cell: str, replicate: int, rank_1: int, rank_2: int, auroc_1: float = 1.0, auroc_2: float = 1.0) -> dict:
    return {
        "cell_id": cell,
        "phase": "power_primary",
        "mode": "orientation",
        "effect": 1.0,
        "replicate": replicate,
        "rank_floor_1": rank_1,
        "rank_floor_2": rank_2,
        "mean_cv_auroc_floor_1": auroc_1,
        "mean_cv_auroc_floor_2": auroc_2,
    }


def test_summary_counts_one_mismatch() -> None:
    rows = [_row("a", r, 3, 3) for r in range(4)] + [_row("b", 0, 2, 2), _row("b", 1, 1, 2, 0.9, 0.95)]
    summary = check.summarize(rows)
    assert summary["count"] == 6
    assert summary["cells"] == 2
    assert summary["mismatches"] == 1
    assert summary["mismatched"] == [{"cell_id": "b", "replicate": 1}]
    assert summary["auroc_differences"] == 1
    assert summary["ranks_floor_1"] == [1, 2, 3]
    assert summary["ranks_floor_2"] == [2, 3]
    assert summary["verdict"] == "ranks_differ"
    assert 1 / 6 < summary["upper_bound_95"] < 1.0


def test_summary_without_mismatch_and_its_bound() -> None:
    rows = [_row(f"c{c}", r, 2, 2) for c in range(19) for r in range(8)]
    summary = check.summarize(rows)
    assert summary["count"] == 152
    assert summary["mismatches"] == 0
    assert summary["verdict"] == "ranks_agree"
    assert summary["upper_bound_95"] == pytest.approx(1 - 0.05 ** (1 / 152))
    assert summary["upper_bound_95"] < 3 / 152


def test_equal_rank_with_differing_auroc_is_not_a_mismatch() -> None:
    summary = check.summarize([_row("a", 0, 3, 3, 0.98, 0.99)])
    assert summary["mismatches"] == 0
    assert summary["auroc_differences"] == 1
    assert summary["verdict"] == "ranks_agree"


def test_duplicate_pairs_and_empty_input_are_refused() -> None:
    with pytest.raises(check.EquivalenceCheckError, match="more than once"):
        check.summarize([_row("a", 0, 2, 2), _row("a", 0, 2, 2)])
    with pytest.raises(check.EquivalenceCheckError, match="no datasets"):
        check.summarize([])


def test_jsonl_conversion_writes_the_csv_and_summary(tmp_path: Path) -> None:
    source = tmp_path / "probe.jsonl"
    records = [
        {"cell_id": "x", "phase": "power_primary", "mode": "shape", "effect": 0.5, "rep": 0, "lv1": 3, "auroc1": 1.0,
         "sec1": 1.0, "lv2": 3, "auroc2": 1.0, "sec2": 1.0},
        {"cell_id": "x", "phase": "power_primary", "mode": "shape", "effect": 0.5, "rep": 1, "lv1": 1, "auroc1": 0.8,
         "sec1": 1.0, "lv2": 2, "auroc2": 0.85, "sec2": 1.0},
    ]
    source.write_text("".join(json.dumps(r) + "\n" for r in records), encoding="utf-8")
    assert check.main(["--from-jsonl", str(source), "--out-dir", str(tmp_path / "out")]) == 0
    with open(tmp_path / "out" / check.CSV_NAME, newline="", encoding="utf-8") as handle:
        table = list(csv.DictReader(handle))
    assert [list(row) for row in table][0] == list(check.CSV_COLUMNS)
    assert [(row["replicate"], row["rank_floor_1"], row["rank_floor_2"]) for row in table] == [
        ("0", "3", "3"),
        ("1", "1", "2"),
    ]
    summary = json.loads((tmp_path / "out" / check.SUMMARY_NAME).read_text(encoding="utf-8"))
    assert summary["mismatches"] == 1 and summary["count"] == 2
    assert summary["provenance"] == {"from_jsonl": str(source)}


def test_jsonl_line_without_a_key_is_refused(tmp_path: Path) -> None:
    source = tmp_path / "probe.jsonl"
    source.write_text(json.dumps({"cell_id": "x", "rep": 0}) + "\n", encoding="utf-8")
    assert check.main(["--from-jsonl", str(source), "--out-dir", str(tmp_path / "out")]) == 1
    assert not (tmp_path / "out").exists()


def test_tasks_follow_the_probe_order() -> None:
    class Cell:
        def __init__(self, cell_id: str, phase: str, mode: str) -> None:
            self.cell_id, self.phase, self.metadata = cell_id, phase, {"trajectory_mode": mode}

    cells = {c.cell_id: c for c in (Cell("t", "type_i_baseline", "none"), Cell("p", "power_primary", "shape"),
                                    Cell("q", "power_primary", "magnitude"))}
    tasks = check.ordered_tasks(cells, [0, 1])
    assert [(cell.cell_id, r) for cell, r in tasks] == [("q", 0), ("q", 1), ("p", 0), ("p", 1), ("t", 0), ("t", 1)]
