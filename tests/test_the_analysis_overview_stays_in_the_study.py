"""The Overview's analysis results are read from the study's own records."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from fastmdxplora.statistics import MEAN_RULES


def _write(root: Path, relative: str, value: object) -> Path:
    path = root / relative
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value), encoding="utf-8")
    return path


def _link(path: Path, target: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    try:
        path.symlink_to(target)
    except OSError:
        pytest.skip("this system cannot create symbolic links")


def _overview(runtime) -> tuple[int, bytes]:
    from tests.test_a_public_dashboard_answers_only_what_it_lists import _ask, _serving

    with _serving(runtime, allow_control=True) as url:
        return _ask(url + "/api/analysis-overview")


def _study(tmp_path: Path) -> Path:
    root = tmp_path / "study"
    _write(root, "manifest.json", {"phases": [{"name": "analysis", "status": "ok"}]})
    _write(root, "analysis/analysis_manifest.json", {"results": {
        "rmsd": {"status": "ok"}, "rg": {"status": "ok"}}})
    _write(root, "analysis/rg/options.json", {"findings": {"mean": {
        "mean": 1.4, "standard_error": 0.1, "n_frames": 200, "effective_samples": 100,
        "rules": MEAN_RULES}}})
    return root


@pytest.mark.parametrize("hidden", [False, True])
def test_a_record_linked_out_of_the_study_or_to_a_hidden_file_is_not_shown(tmp_path, hidden):
    from fastmdxplora.gui.exploration import DashboardRuntime

    root = _study(tmp_path)
    secret = {"findings": {"mean": {"mean": 9.9, "n_frames": 200,
                                    "not_a_measurement": "PRIVATE_OUTSIDE_TEXT"}}}
    target = (_write(root, ".private/options.json", secret) if hidden
              else _write(tmp_path, "outside/options.json", secret))
    _link(root / "analysis/rmsd/options.json", target)
    runtime = DashboardRuntime(workspace_root=root, exploration_root=tmp_path, active_root=root)
    status, body = _overview(runtime)
    assert status == 200
    assert b"PRIVATE_OUTSIDE_TEXT" not in body
    rows = {row["analysis"]: row for row in json.loads(body)["rows"]}
    # The study's own record is still read.
    assert rows["rg"]["mean"]["determined"] is True


def test_a_link_inside_the_study_is_still_followed(tmp_path):
    from fastmdxplora.gui.exploration import DashboardRuntime

    root = _study(tmp_path)
    kept = _write(root, "analysis/rmsd/kept.json", {"findings": {"mean": {
        "mean": 0.2, "n_frames": 200, "not_a_measurement": "Recorded inside the study."}}})
    _link(root / "analysis/rmsd/options.json", kept)
    runtime = DashboardRuntime(workspace_root=root, exploration_root=tmp_path, active_root=root)
    status, body = _overview(runtime)
    assert status == 200
    assert b"Recorded inside the study." in body


def test_no_study_open_answers_an_empty_overview(tmp_path):
    from fastmdxplora.gui.exploration import DashboardRuntime

    runtime = DashboardRuntime(workspace_root=tmp_path, exploration_root=tmp_path,
                               active_root=None)
    status, body = _overview(runtime)
    assert status == 200
    assert json.loads(body)["rows"] == []
