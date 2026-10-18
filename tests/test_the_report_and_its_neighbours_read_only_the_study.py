"""The Report page, the methods, the colours by residue, the charts and
the run's narration read only the study's own records.

These routes read through a link out of the study: the Report page showed
an outside report.md and its figure panel an outside options.json, and the
methods, the per-residue values, the charts and the events read what such
links led to. Each now reads in a held read (`server._confined`), names the
records it did not read, and says so; a list read across the study is not
given in part.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from tests.test_the_study_pages_read_only_the_study import (
    _OPENED,
    _WATCHED,
    _a_study_carried_on,
    _ask,
    _every_record_linked_out,
    _link,
    _watch,
    _write,
)

PUBLIC = pytest.mark.parametrize("public", [False, True], ids=["loopback", "public"])
ROUTES = ("/api/report", "/api/methods", "/api/residue-values", "/api/metrics", "/api/events")
MARK = "OUTSIDE-REPORT-not-for-any-page"


def _with_its_report(tmp_path: Path) -> Path:
    """A study carried on, with what each of these routes reads: its
    report, a per-residue table, and each piece's charts and narration."""
    root = _a_study_carried_on(tmp_path)
    _write(root, "report/report.md", "# The study\n\n_Generated: today_\n\nWhat it found.\n")
    _write(root, "analysis/rmsf/options.json", {"analysis": "rmsf", "options": {}})
    _write(root, "analysis/rmsf/rmsf.dat", "# residue rmsf\n1 0.10\n2 0.12\n3 0.20\n")
    for piece in ("", "segment-001/"):
        # Counted, so the pieces' samples are read as one run's.
        _write(root, f"{piece}simulation/live_status.json", {
            "status": "completed", "stage": "production", "current_step": 10000,
            "total_planned_steps": 10000, "current_frame_count": 100, "timestep_fs": 2.0,
            "last_update_timestamp": f"2026-10-0{2 if piece else 1}T00:00:00+00:00"})
        _write(root, f"{piece}simulation/live_metrics.csv",
               "stage,simulation_time_ns,temperature\nproduction,0.001,300.0\n")
        _write(root, f"{piece}simulation/live_events.log",
               "2026-10-01T00:00:00+00:00\tinfo\tProduction began.\n")
    return root


@PUBLIC
@pytest.mark.parametrize("route", ROUTES)
def test_no_record_is_read_through_a_link_out(tmp_path, route, public):
    """Every record of the study linked out of it: none is opened where it
    leads, and the answer names what was not read."""
    root = _with_its_report(tmp_path)
    outside = _every_record_linked_out(root, tmp_path)
    _watch(outside)
    try:
        status, body = _ask(root, tmp_path, route, public=public)
        opened = list(_OPENED)
    finally:
        _WATCHED.clear()
    assert status == 200
    assert opened == []
    assert json.loads(body)["unavailable"]


@PUBLIC
def test_a_report_linked_out_is_said_not_read(tmp_path, public):
    root = _with_its_report(tmp_path)
    target = tmp_path / "outside" / "report.md"
    target.parent.mkdir(parents=True)
    target.write_text(f"# {MARK}\n", encoding="utf-8")
    (root / "report" / "report.md").unlink()
    _link(root / "report" / "report.md", target)

    status, body = _ask(root, tmp_path, "/api/report", public=public)
    assert MARK not in body.decode("utf-8")
    answer = json.loads(body)
    assert answer["ok"] is False and answer["reason"] != "no report yet"
    assert answer["not_read"].startswith("report/report.md was not read: it is a link to a "
                                         "file outside the study")


@PUBLIC
@pytest.mark.parametrize("route,key", [("/api/metrics", "metrics"), ("/api/events", "events")])
def test_a_piece_s_record_not_read_gives_none_of_the_rest(tmp_path, route, key, public):
    """A piece's charts or narration linked out: not the first piece's
    alone, said as the whole; none, and why."""
    root = _with_its_report(tmp_path)
    name = "live_metrics.csv" if key == "metrics" else "live_events.log"
    _every_record_linked_out(root, tmp_path, named=name)
    (root / "simulation" / name).unlink()
    _write(root, f"simulation/{name}", (tmp_path / "outside" / "simulation" / name)
           .read_text(encoding="utf-8"))

    answer = json.loads(_ask(root, tmp_path, route, public=public)[1])
    assert answer[key] == []
    assert answer["unavailable"] == [f"segment-001/simulation/{name}"]
    assert answer["not_read"].startswith(f"segment-001/simulation/{name} was not read")


@PUBLIC
def test_every_record_read_the_answers_are_as_before(tmp_path, public):
    """No link: each route answers as it did, with nothing said not read."""
    root = _with_its_report(tmp_path)
    for route in ROUTES:
        answer = json.loads(_ask(root, tmp_path, route, public=public)[1])
        assert "unavailable" not in answer and "not_read" not in answer, route
    assert json.loads(_ask(root, tmp_path, "/api/report", public=public)[1])["ok"] is True
    values = json.loads(_ask(root, tmp_path, "/api/residue-values", public=public)[1])
    assert [entry["key"] for entry in values["properties"]] == ["rmsf"]
    assert len(json.loads(_ask(root, tmp_path, "/api/events", public=public)[1])["events"]) == 2
