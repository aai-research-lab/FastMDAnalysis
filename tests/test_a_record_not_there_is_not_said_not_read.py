"""A record that is not there is not said not read.

A study of several runs, one of them carried on, was said on the Report page
and in the colours by residue to have a record "not read": the process
record of the piece being run, which a held read asks for whether or not it
is there, and which was admitted only in a study of one run. Such a study
has no link, and answered as if it had one. A run's own process record is
admitted in a study of several as in a study of one, and a record that is
not there is missing, not refused: only a link or a hidden file that is
there is said not read.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from fastmdxplora.orchestrator import RUN_PROCESS_FILE
from tests.test_the_study_pages_read_only_the_study import _ask, _link, _write

PUBLIC = pytest.mark.parametrize("public", [False, True], ids=["loopback", "public"])
ROUTES = ("/api/report", "/api/residue-values")


def _runs(tmp_path: Path) -> Path:
    """A study of two runs, the first carried on in `segment-001`, whose
    record says it is going in production there."""
    from datetime import datetime, timezone

    root = tmp_path / "study"
    _write(root, "batch_manifest.json", {
        "planned": [{"run_id": "r1"}, {"run_id": "r2"}],
        "runs": [{"run_id": "r1", "status": "ok"}, {"run_id": "r2", "status": "ok"}]})
    now = datetime.now(timezone.utc).isoformat()
    for run in ("r1", "r2"):
        _write(root, f"runs/{run}/simulation/live_status.json", {
            "status": "completed", "stage": "production", "current_step": 1000,
            "total_planned_steps": 1000, "last_update_timestamp": "2026-10-01T00:00:00+00:00"})
    _write(root, "runs/r1/segment-001/simulation/live_status.json", {
        "status": "running", "stage": "production", "current_step": 1500,
        "total_planned_steps": 2000, "last_update_timestamp": now})
    return root


def _unheld(root: Path, route: str) -> dict:
    """What the route's reader answers outside a held read."""
    from fastmdxplora.gui.by_residue import values_by_residue
    from fastmdxplora.gui.report_page import report_payload

    read = report_payload if route == "/api/report" else values_by_residue
    return json.loads(json.dumps(read(root)))


@PUBLIC
@pytest.mark.parametrize("route", ROUTES)
def test_a_process_record_not_there_is_not_said_not_read(tmp_path, route, public):
    root = _runs(tmp_path)
    answer = json.loads(_ask(root, tmp_path, route, public=public)[1])
    assert "unavailable" not in answer and "not_read" not in answer
    assert answer == _unheld(root, route)


@PUBLIC
def test_a_run_s_own_process_record_is_read(tmp_path, public):
    """The piece's process record there, its process gone: the run is read
    as interrupted, as outside a held read, and nothing is said not read."""
    from fastmdxplora.orchestrator import this_machine

    root = _runs(tmp_path)
    here = this_machine()
    record = {**here, "pid": 4_194_000, "started_at": "2026-10-01T00:00:00+00:00",
              "argv": ["fastmdx", "resume"]}
    if here.get("boot"):
        record["boot"] = "a-boot-before-this-one"
    _write(root, f"runs/r1/segment-001/{RUN_PROCESS_FILE}", record)
    answer = json.loads(_ask(root, tmp_path, "/api/report", public=public)[1])
    assert "unavailable" not in answer and "not_read" not in answer
    assert answer == _unheld(root, "/api/report")


def test_a_run_s_own_process_record_says_how_it_ended(tmp_path):
    """Read in a held read as outside one: the piece's process gone, the
    run is interrupted, not running."""
    from fastmdxplora.gui.exploration import runs_of_a_study
    from fastmdxplora.gui.server import _confined
    from fastmdxplora.orchestrator import this_machine

    root = _runs(tmp_path)
    here = this_machine()
    record = {**here, "pid": 4_194_000, "started_at": "2026-10-01T00:00:00+00:00",
              "argv": ["fastmdx", "resume"]}
    if here.get("boot"):
        record["boot"] = "a-boot-before-this-one"
    _write(root, f"runs/r1/segment-001/{RUN_PROCESS_FILE}", record)
    held = _confined(root, lambda inside: {"runs": runs_of_a_study(inside)})
    assert "unavailable" not in held
    assert {run["run_id"]: run["state"] for run in held["runs"]}["r1"] == "interrupted"
    assert held["runs"] == runs_of_a_study(root)


def test_a_result_linked_to_a_run_s_process_record_is_refused(tmp_path):
    """Admitted by its own name only: a result of a run linked to the run's
    process record is refused, not read as that result."""
    from fastmdxplora import record_access
    from fastmdxplora.gui.server import _confined

    root = _runs(tmp_path)
    _write(root, f"runs/r1/{RUN_PROCESS_FILE}", {"pid": 1})
    _link(root / "runs/r1/analysis/rmsd.dat", root / "runs/r1" / RUN_PROCESS_FILE)

    def read(inside: Path) -> dict:
        try:
            return {"said": record_access.read_text(inside / "runs/r1/analysis/rmsd.dat")}
        except OSError:
            return {"said": None}

    answer = _confined(root, read)
    assert answer["said"] is None and answer["unavailable"] == ["runs/r1/analysis/rmsd.dat"]
    own = _confined(root, lambda inside: {"said": record_access.read_text(
        inside / "runs/r1" / RUN_PROCESS_FILE)})
    assert json.loads(own["said"]) == {"pid": 1} and "unavailable" not in own


def test_only_a_record_that_is_there_is_refused(tmp_path):
    """A hidden file that is not there is missing; one that is there, and a
    link out (even to nothing), are refused and named."""
    from fastmdxplora import record_access

    root = tmp_path / "study"
    _write(root, ".hidden-there", "secret")
    _link(root / "leads-out", tmp_path / "nothing-there")

    def admitted(path: Path) -> bool:
        return not any(part.startswith(".") for part in path.relative_to(root).parts)

    with record_access.confined_records(root, admitted) as scope:
        with pytest.raises(FileNotFoundError):
            record_access.read_text(root / ".hidden-not-there")
        assert record_access.can_read(root / ".hidden-not-there") is False
        assert record_access.is_file(root / "segment-004" / ".hidden-not-there") is False
        assert scope.denied == set()
        with pytest.raises(PermissionError):
            record_access.read_text(root / ".hidden-there")
        with pytest.raises(PermissionError):
            record_access.read_text(root / "leads-out")
    assert scope.denied == {".hidden-there", "leads-out"}
