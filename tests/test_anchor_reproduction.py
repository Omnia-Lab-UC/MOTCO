"""``scripts/anchor_reproduction.py`` on synthetic merged JSONL.

The by-column mode pairs zero-effect anchors on (design column, replicate
index), the column read from ``cell_metadata.design_point`` (``power_primary``
is the baseline column, ``type_i_baseline`` anchors are excluded), and reports
one summary line per column. Without the flag the single-anchor path writes
exactly the columns and summary it wrote before the mode existed.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
try:
    import anchor_reproduction as ar
finally:
    sys.path.pop(0)

N = "generator.n_samples"
LAYERS = "evaluation.integration_params.layers"
TWO = ["methylation", "expression"]
BASELINE = {N: 200, LAYERS: TWO}
COLUMNS = [{N: 80, LAYERS: TWO}, BASELINE, {N: 200, LAYERS: None}]

#: The single-anchor CSV columns, exactly as written before the by-column mode existed.
PRE_CHANGE_CSV_COLUMNS = (
    "cell_id_candidate,replicate_index,status_candidate,generator_seed_candidate,selected_lv_candidate,"
    "permutations_candidate,observed_delta_candidate,p_delta_candidate,observed_angle_candidate,p_angle_candidate,"
    "observed_shape_candidate,p_shape_candidate,cell_id_reference,status_reference,generator_seed_reference,"
    "selected_lv_reference,permutations_reference,observed_delta_reference,p_delta_reference,"
    "observed_angle_reference,p_angle_reference,observed_shape_reference,p_shape_reference,seed_match,lv_match,"
    "abs_diff_delta,delta_match,p_delta_match,abs_diff_angle,angle_match,p_angle_match,abs_diff_shape,shape_match,"
    "p_shape_match,p_values_comparable,all_match"
).split(",")


def _anchor(
    cell_id: str,
    phase: str,
    index: int,
    *,
    permutations: int,
    point: dict | None,
    seed_offset: int = 0,
    p_offset: float = 0.0,
) -> dict:
    """A deterministic anchor record: seeds and statistics are functions of the replicate index only."""

    record = {
        "cell_id": cell_id,
        "phase": phase,
        "replicate_index": index,
        "replicate_seed": 1000 + index,
        "generator_seed": 5000 + index + seed_offset,
        "evaluation_seed": 3,
        "parameter_signature": f"sig-{cell_id}",
        "status": "completed",
        "p_values": {
            s: round(((index * 7 + k * 3) % (permutations + 1) + 1) / (permutations + 1) + p_offset, 6)
            for k, s in enumerate(("delta", "angle", "shape"))
        },
        "pair_statistics": {"delta": 0.1 * index, "angle": 30.0 + index, "shape": 0.01 * index},
        "runtime_metadata": {"permutations": permutations, "n_jobs": 1},
        "cell_metadata": {"trajectory_mode": "none", "effect_size": 0.0, "zero_effect_anchor": True},
        "integration_metadata": {"integration_method": "pls", "selected_lv": 2 + (index % 2)},
    }
    if point is not None:
        record["cell_metadata"]["design_point"] = point
    return record


def _power_cell(cell_id: str, phase: str, index: int, point: dict | None) -> dict:
    record = _anchor(cell_id, phase, index, permutations=99, point=point)
    record["cell_metadata"].update({"trajectory_mode": "magnitude", "effect_size": 0.5, "zero_effect_anchor": False})
    return record


def _column_anchor_id(prefix: str, point: dict) -> str:
    blocks = 3 if point[LAYERS] is None else len(point[LAYERS])
    return f"{prefix}-anchor-{point[N]}-{blocks}"


def _design_grid_run(
    prefix: str, *, n_replicates: int, permutations: int, columns: list[dict] = COLUMNS, type_i_seed_offset: int = 0
) -> list[dict]:
    records = []
    for point in columns:
        phase = "power_primary" if point == BASELINE else "power_design"
        cell_id = _column_anchor_id(prefix, point)
        records.extend(_anchor(cell_id, phase, i, permutations=permutations, point=point) for i in range(n_replicates))
        records.extend(_power_cell(f"{prefix}-magnitude-{point[N]}", phase, i, point) for i in range(n_replicates))
    # The Type I baseline anchor: a different draw, never paired.
    records.extend(
        _anchor(f"{prefix}-type-i-anchor", "type_i_baseline", i, permutations=permutations, point=None,
                seed_offset=type_i_seed_offset)
        for i in range(n_replicates)
    )
    return records


def _write(path: Path, records: list[dict]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for record in records:
            handle.write(json.dumps(record, sort_keys=True) + "\n")
    return path


def _run(capsys, *argv: str) -> tuple[int, list[str]]:
    code = ar.main(list(argv))
    return code, capsys.readouterr().out.splitlines()


def _by_column(capsys, candidate: Path, reference: Path, out: Path) -> tuple[int, list[str]]:
    return _run(capsys, "--by-column", "--candidate", str(candidate), "--reference", str(reference), "--out", str(out))


def _column_lines(lines: list[str]) -> dict[str, str]:
    """Per-column summary lines keyed by the column's JSON design point (the key itself holds colons)."""

    out = {}
    for line in lines:
        if line.startswith("  column "):
            body = line[len("  column ") :]
            end = body.rindex("}") + 1
            out[body[:end]] = body[end:].removeprefix(" [baseline]").removeprefix(": ")
    return out


def test_sibling_runs_match_at_every_column(tmp_path: Path, capsys) -> None:
    candidate = _write(tmp_path / "magnitude.jsonl", _design_grid_run("mag", n_replicates=5, permutations=999))
    reference = _write(
        tmp_path / "geometry.jsonl",
        _design_grid_run("geo", n_replicates=5, permutations=999, type_i_seed_offset=777),
    )
    out = tmp_path / "siblings.csv"
    code, lines = _by_column(capsys, candidate, reference, out)
    assert code == 0
    frame = pd.read_csv(out)
    assert len(frame) == 3 * 5
    assert frame["all_match"].all() and frame["p_values_comparable"].all()
    assert set(frame["design_column"]) == {ar.column_key(point) for point in COLUMNS}
    # The column coordinates are written, one CSV column per axis, ahead of the replicate index.
    assert list(frame.columns[:4]) == ["design_column", LAYERS, N, "replicate_index"]
    assert set(frame[N]) == {"80", "200"} or set(frame[N]) == {80, 200}
    # The Type I baseline anchors (different seeds on each side) are not paired.
    assert "type-i" not in " ".join(frame["cell_id_candidate"].astype(str))
    assert lines[0].startswith(
        "candidate 15 anchor records over 3 design column(s) vs reference 15 anchor records over 3"
    )
    per_column = _column_lines(lines)
    assert set(per_column) == {ar.column_key(point) for point in COLUMNS}
    for rest in per_column.values():
        assert rest.startswith("5 pairs (candidate 5, reference 5); seeds identical 5/5; rank identical 5/5")
        assert "p-values identical 5/5" in rest and rest.endswith("all comparable fields match 5/5")
    assert sum(" [baseline]" in line for line in lines) == 1
    assert lines[1].startswith(f"  column {ar.column_key(BASELINE)} [baseline]")
    assert "  all comparable fields match: 15/15 (paired columns only)" in lines


def test_pilot_run_matches_seeds_and_statistics_with_p_values_not_comparable(tmp_path: Path, capsys) -> None:
    candidate = _write(tmp_path / "slice.jsonl", _design_grid_run("slice", n_replicates=6, permutations=999))
    reference = _write(tmp_path / "pilot.jsonl", _design_grid_run("pilot", n_replicates=3, permutations=199))
    out = tmp_path / "pilot.csv"
    code, lines = _by_column(capsys, candidate, reference, out)
    assert code == 0
    frame = pd.read_csv(out)
    assert len(frame) == 3 * 3  # replicate indices 0-2 only
    assert set(frame["replicate_index"]) == {0, 1, 2}
    assert frame["seed_match"].all() and frame["delta_match"].all() and frame["angle_match"].all()
    assert not frame["p_values_comparable"].any()
    assert not frame["p_delta_match"].all()  # the p-values do differ at 999 vs 199 permutations ...
    assert frame["all_match"].all()  # ... and are excluded from the verdict
    for rest in _column_lines(lines).values():
        assert "p-values not comparable ([999] vs [199])" in rest
        assert rest.endswith("all comparable fields match 3/3")
    assert "  p-values not comparable on some pairs (permutation counts differ); excluded from the verdict" in lines


def test_column_missing_on_one_side_is_listed_with_zero_pairs(tmp_path: Path, capsys) -> None:
    candidate = _write(tmp_path / "candidate.jsonl", _design_grid_run("cand", n_replicates=4, permutations=99))
    reference = _write(
        tmp_path / "reference.jsonl", _design_grid_run("ref", n_replicates=4, permutations=99, columns=COLUMNS[:2])
    )
    out = tmp_path / "missing.csv"
    code, lines = _by_column(capsys, candidate, reference, out)
    assert code == 0  # the verdict covers the paired columns only
    frame = pd.read_csv(out)
    assert len(frame) == 2 * 4
    per_column = _column_lines(lines)
    assert per_column[ar.column_key(COLUMNS[2])] == "0 pairs (candidate 4 records, reference 0 records)"
    assert "  all comparable fields match: 8/8 (paired columns only)" in lines
    assert lines[0].endswith("8 (design column, replicate index) pairs over 2 paired column(s)")


def test_by_column_mismatch_fails_the_verdict(tmp_path: Path, capsys) -> None:
    candidate = _write(tmp_path / "candidate.jsonl", _design_grid_run("cand", n_replicates=3, permutations=99))
    altered = _design_grid_run("ref", n_replicates=3, permutations=99)
    for record in altered:
        if record["cell_metadata"].get("zero_effect_anchor") and record["phase"] == "power_design":
            record["generator_seed"] += 1
    reference = _write(tmp_path / "reference.jsonl", altered)
    out = tmp_path / "altered.csv"
    code, lines = _by_column(capsys, candidate, reference, out)
    assert code == 1
    frame = pd.read_csv(out)
    assert int(frame["seed_match"].sum()) == 3  # only the baseline column's seeds still match
    assert _column_lines(lines)[ar.column_key(BASELINE)].endswith("all comparable fields match 3/3")
    assert _column_lines(lines)[ar.column_key(COLUMNS[0])].endswith("all comparable fields match 0/3")


def test_by_column_refuses_cell_selection_and_crowded_columns(tmp_path: Path) -> None:
    run = _design_grid_run("cand", n_replicates=2, permutations=99)
    candidate = _write(tmp_path / "candidate.jsonl", run)
    with pytest.raises(SystemExit):
        ar.main(["--by-column", "--candidate", str(candidate), "--reference", str(candidate),
                 "--candidate-cell", "x", "--out", str(tmp_path / "o.csv")])
    crowded = run + [_anchor("cand-second-anchor", "power_design", 0, permutations=99, point=COLUMNS[0])]
    crowded_path = _write(tmp_path / "crowded.jsonl", crowded)
    with pytest.raises(SystemExit, match="several anchor cells at one design column"):
        ar.main(["--by-column", "--candidate", str(crowded_path), "--reference", str(candidate),
                 "--out", str(tmp_path / "o.csv")])


def test_default_mode_writes_the_pre_change_csv_and_summary(tmp_path: Path, capsys) -> None:
    candidate = _write(
        tmp_path / "candidate.jsonl",
        [_anchor("cand-anchor", "power_primary", i, permutations=99, point=None) for i in range(6)],
    )
    reference_records = [_anchor("ref-anchor", "power_primary", i, permutations=99, point=None) for i in range(6)]
    reference_records.append(_power_cell("ref-power", "power_primary", 0, None))
    reference = _write(tmp_path / "reference.jsonl", reference_records)
    out = tmp_path / "default.csv"
    code, lines = _run(capsys, "--candidate", str(candidate), "--reference", str(reference), "--out", str(out))
    assert code == 0
    assert list(pd.read_csv(out).columns) == PRE_CHANGE_CSV_COLUMNS
    assert lines == [
        "candidate cand-anchor (6 records) vs reference ref-anchor (6 records); 6 replicate indices in common",
        "  seeds identical: 6/6; selected rank identical: 6/6",
        "  delta: observed within rtol 1e-06 on 6/6 (max abs diff 0); p-values identical on 6/6",
        "  angle: observed within rtol 1e-06 on 6/6 (max abs diff 0); p-values identical on 6/6",
        "  shape: observed within rtol 1e-06 on 6/6 (max abs diff 0); p-values identical on 6/6",
        "  all comparable fields match: 6/6",
        f"Wrote {out}",
    ]
