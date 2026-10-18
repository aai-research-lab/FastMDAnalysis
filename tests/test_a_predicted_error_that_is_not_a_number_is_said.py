"""Replicas whose runs gave no error of their own are not said to agree with one.

A run too short to estimate its own error records NaN. NaN compared with
anything is false, so the spread of the replicas fell through to "within the
factor of 2 ... against a predicted nan ... supported by repeating it".
"""

from __future__ import annotations

import json
import math
from pathlib import Path

from fastmdxplora.batch.aggregate import aggregate_members


def _seed_campaign(tmp_path: Path, errors: list[float | None]) -> Path:
    runs = []
    for i, se in enumerate(errors):
        out = tmp_path / "runs" / f"seed{i}"
        (out / "analysis" / "rg").mkdir(parents=True)
        mean = {"mean": 1.00 + 0.01 * i}
        if se is not None:
            mean["standard_error"] = se
        # Written as the analyses write it: json.dumps writes NaN as NaN.
        (out / "analysis" / "rg" / "options.json").write_text(
            json.dumps({"analysis": "rg", "findings": {"mean": mean}}),
            encoding="utf-8")
        runs.append({"run_id": f"seed{i}", "output_dir": str(out),
                     "status": "ok", "system": "a.pdb",
                     "sweep_values": {"simulation.random_seed": i}})
    (tmp_path / "batch_manifest.json").write_text(json.dumps(
        {"runs": runs, "sweep": {"simulation.random_seed": [0, 1]}}),
        encoding="utf-8")
    return tmp_path


def test_no_finite_prediction_is_said_as_none(tmp_path):
    entry = aggregate_members(
        _seed_campaign(tmp_path, [float("nan"), float("nan")]))["analyses"]["rg"]
    said = entry["calibration"]

    assert "nan" not in said.lower()
    assert "supported" not in said
    assert "within the factor" not in said
    assert "No run gave an error of its own" in said
    assert entry["predicted_standard_error"] is None
    assert entry["ratio_observed_to_predicted"] is None
    assert math.isfinite(entry["observed_standard_error"])
    # The record stays JSON a browser can read: NaN is not JSON.
    json.dumps(entry, allow_nan=False)


def test_a_prediction_from_some_runs_says_how_many(tmp_path):
    entry = aggregate_members(
        _seed_campaign(tmp_path, [0.004, float("nan"), 0.006]))["analyses"]["rg"]
    said = entry["calibration"]

    assert entry["predicted_standard_error"] == 0.005
    assert "nan" not in said.lower()
    assert "from the 2 of 3 runs that gave an error of their own" in said


def test_a_zero_error_is_not_a_prediction(tmp_path):
    entry = aggregate_members(
        _seed_campaign(tmp_path, [0.0, 0.0]))["analyses"]["rg"]

    assert "supported" not in entry["calibration"]
    assert entry["predicted_standard_error"] is None
