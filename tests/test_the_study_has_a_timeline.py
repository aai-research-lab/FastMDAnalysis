"""The study's Timeline: what happened to it, newest first, from its records.

A study's history was in pieces nobody read together: the Manifest's record
of each phase as it last ran, the record and output of a phase a run again
replaced (`previous/<phase>`) and the pieces an extended study was carried
on in. The side panel's Timeline lists them as an editor lists a file's
history: when, how long, which version on which computer, and what each
changed where the records say (the means an analysis run again moved, the
settings a piece was carried on with). Who asked for a step is not
recorded, and is not said.
"""

from __future__ import annotations

import json
import shutil
import urllib.request
from pathlib import Path

import pytest

from tests.test_a_study_of_runs_is_shown_as_one import _recorded
from tests.test_the_frame_holds_the_study import _a_finished_study

MADE = {"version": "2.6.1", "host": "aailab01"}


def _phase(name: str, start: str, end: str) -> dict:
    return {"name": name, "status": "ok", "started_at": f"2026-10-18T{start}+00:00",
            "finished_at": f"2026-10-18T{end}+00:00", "produced_by": MADE}


def _options(folder: Path, name: str, mean: dict) -> None:
    (folder / name).mkdir(parents=True, exist_ok=True)
    (folder / name / "options.json").write_text(json.dumps(
        {"analysis": name, "options": {}, "findings": {"mean": mean}}), encoding="utf-8")


@pytest.fixture()
def study(tmp_path) -> Path:
    root = tmp_path / "study"
    shutil.copytree(_a_finished_study(), root)
    (root / "manifest.json").write_text(json.dumps({"phases": [
        _phase("setup", "07:00:00", "07:00:04"), _phase("simulation", "07:00:04", "07:30:04"),
        _phase("analysis", "09:00:00", "09:00:12")]}), encoding="utf-8")
    (root / "analysis" / "analysis_manifest.json").write_text(
        json.dumps({"plan": ["rmsd", "rg"]}), encoding="utf-8")
    _options(root / "analysis", "rmsd", _recorded(0.16))
    _options(root / "analysis", "rg", _recorded(1.2))
    before = root / "previous" / "analysis"
    _options(before, "rmsd", _recorded(0.15))
    _options(before, "rg", _recorded(1.2))
    _options(before, "hbonds", _recorded(12.0))
    (before / "phase_record.json").write_text(json.dumps(
        _phase("analysis", "07:30:04", "07:31:00")), encoding="utf-8")
    report = root / "previous" / "report"
    report.mkdir(parents=True)
    (report / "phase_record.json").write_text(json.dumps(
        _phase("report", "07:31:00", "07:31:30")), encoding="utf-8")
    piece = root / "segment-001"
    (piece / "simulation").mkdir(parents=True)
    (piece / "manifest.json").write_text(json.dumps({"phases": [
        _phase("simulation", "08:00:00", "08:20:00")]}), encoding="utf-8")
    config = (root / "resolved_config.yml").read_text(encoding="utf-8")
    (piece / "resolved_config.yml").write_text(
        config.replace("duration_ns: 1.0", "duration_ns: 1.0\n  extra_ns: 2.0"), encoding="utf-8")
    return root


def test_the_steps_are_listed_newest_first_with_what_each_changed(study):
    from fastmdxplora.gui.study_timeline import timeline_of

    events = timeline_of(study)["events"]
    assert [e["said"] for e in events] == [
        "Analysed", "Carried on, piece 1",
        "Report written, then set aside: written from what was run again",
        "Analysed, then replaced", "Simulated", "Set up"]
    analysed = events[0]
    assert analysed["again"] is True and analysed["seconds"] == 12.0
    assert (analysed["version"], analysed["host"]) == ("2.6.1", "aailab01")
    moved = analysed["changes"]["moved"]
    assert [m["analysis"] for m in moved] == ["rmsd"]
    assert moved[0]["before"]["mean"] == pytest.approx(0.15, abs=1e-3)
    assert moved[0]["after"]["mean"] == pytest.approx(0.16, abs=1e-3)
    assert analysed["changes"]["kept"] == 1 and analysed["changes"]["dropped"] == ["hbonds"]
    assert analysed["detail"] == "2 analyses: rmsd, rg"
    piece = events[1]
    assert piece["seconds"] == 1200.0
    assert {d["setting"] for d in piece["settings"]} >= {"simulation.extra_ns"}
    assert events[2]["where"] == "previous/report"
    assert all("person" not in e["said"].lower() and "agent" not in e["said"].lower()
               for e in events)


def test_a_folder_that_is_no_study_is_said(tmp_path):
    from fastmdxplora.gui.study_timeline import timeline_of

    assert timeline_of(tmp_path / "nothing")["ok"] is False
    assert timeline_of(tmp_path) == {"ok": True, "events": []}


def test_the_route_is_answered_here_only(study):
    from fastmdxplora.gui.server import GETS_ANSWERED_BEYOND_LOOPBACK, start_dashboard_session

    assert "/api/study-timeline" not in GETS_ANSWERED_BEYOND_LOOPBACK
    session = start_dashboard_session(output=str(study), host="127.0.0.1", port=0)
    try:
        said = json.loads(urllib.request.urlopen(
            session.url.rstrip("/") + "/api/study-timeline", timeout=30).read())
    finally:
        session.server.shutdown()
    assert said["ok"] is True and len(said["events"]) == 6


def test_the_side_panel_lists_them(study):
    pytest.importorskip("playwright.sync_api")
    from playwright.sync_api import sync_playwright

    from fastmdxplora.gui.server import start_dashboard_session

    session = start_dashboard_session(output=str(study), host="127.0.0.1", port=0)
    try:
        with sync_playwright() as pw:
            browser = pw.chromium.launch(args=["--enable-unsafe-swiftshader"])
            page = browser.new_page(viewport={"width": 1600, "height": 900})
            page.set_default_timeout(60000)
            errors: list[str] = []
            page.on("pageerror", lambda error: errors.append(str(error)))
            page.goto(session.url + "#overview", wait_until="domcontentloaded")
            page.wait_for_selector("body.state-ready")
            page.evaluate("() => window.FastMDXFrame.setCollapsed(false)")
            page.click("#timeline-tab")
            page.wait_for_selector("#study-timeline-list .tl-event")
            page.click("#study-timeline-list .tl-changes > summary")
            said = page.evaluate("""() => ({
                steps: [...document.querySelectorAll('.tl-event .tl-said')].map((n) => n.textContent),
                changes: document.querySelector('.tl-changes > summary').textContent,
                moved: [...document.querySelectorAll('.tl-table tr[data-analysis]')]
                    .map((r) => r.dataset.analysis),
                set_aside: document.querySelectorAll('.tl-event[data-replaced="true"]').length,
            })""")
            browser.close()
    finally:
        session.server.shutdown()
    assert errors == []
    assert said["steps"][0] == "Analysed again" and said["steps"][-1] == "Set up"
    assert said["changes"] == "1 mean moved, 1 the same, 1 not run again"
    assert said["moved"] == ["rmsd"] and said["set_aside"] == 2
