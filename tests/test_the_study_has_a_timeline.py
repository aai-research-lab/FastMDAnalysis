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


# ---------------------------------------------------------------------------
# Read held to the study (`server._confined`): a record linked out of it is
# not read, and said so, never as a step or a mean the study does not have.

PUBLIC = pytest.mark.parametrize("public", [False, True], ids=["loopback", "public"])


def _timeline(study: Path, workspace: Path, outside: Path, *, public: bool) -> tuple[int, dict]:
    from tests.test_the_study_pages_read_only_the_study import _ask, _OPENED, _watch, _WATCHED

    _watch(outside)
    try:
        status, body = _ask(study, workspace, "/api/study-timeline", public=public)
        opened = sorted(set(_OPENED))
    finally:
        _WATCHED.clear()
    assert opened == []
    return status, json.loads(body)


def _linked_out(record: Path, outside: Path) -> None:
    from tests.test_the_study_pages_read_only_the_study import _link

    outside.parent.mkdir(parents=True, exist_ok=True)
    record.rename(outside)
    _link(record, outside)


@PUBLIC
def test_a_mean_s_record_linked_out_is_said_not_read_not_added_or_dropped(study, tmp_path,
                                                                           public):
    outside = tmp_path / "elsewhere"
    _linked_out(study / "analysis" / "rmsd" / "options.json", outside / "now.json")
    _linked_out(study / "previous" / "analysis" / "hbonds" / "options.json",
                outside / "before.json")
    status, said = _timeline(study, tmp_path, outside, public=public)
    if public:
        # Not answered beyond loopback at all.
        assert status == 403
        return
    assert status == 200 and said["ok"] is True
    changes = next(e for e in said["events"] if e["said"] == "Analysed")["changes"]
    assert changes["not_read"] == ["hbonds", "rmsd"]
    assert changes["moved"] == [] and changes["added"] == [] and changes["dropped"] == []
    assert said["unavailable"] == ["analysis/rmsd/options.json",
                                   "previous/analysis/hbonds/options.json"]
    assert said["not_read"].startswith("analysis/rmsd/options.json was not read")


def test_a_piece_s_settings_linked_out_are_said_not_read(study, tmp_path):
    outside = tmp_path / "elsewhere"
    _linked_out(study / "segment-001" / "resolved_config.yml", outside / "resolved_config.yml")
    status, said = _timeline(study, tmp_path, outside, public=False)
    piece = next(e for e in said["events"] if e["kind"] == "piece")
    assert status == 200
    assert piece["settings_not_read"] is True and piece["settings"] == []
    assert "segment-001/resolved_config.yml was not read" in said["not_read"]


def test_the_side_panel_says_a_record_not_read(study, tmp_path):
    pytest.importorskip("playwright.sync_api")
    from playwright.sync_api import sync_playwright

    from fastmdxplora.gui.server import start_dashboard_session

    _linked_out(study / "analysis" / "rmsd" / "options.json", tmp_path / "elsewhere" / "now.json")
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
                line: document.getElementById('study-timeline-said').textContent,
                changes: document.querySelector('.tl-changes > summary').textContent,
                moved: document.querySelectorAll('.tl-table tr[data-analysis]').length,
            })""")
            browser.close()
    finally:
        session.server.shutdown()
    assert errors == []
    assert said["line"].startswith("analysis/rmsd/options.json was not read")
    assert said["changes"] == "1 the same, 1 not run again, 1 not read"
    assert said["moved"] == 0


# ---------------------------------------------------------------------------
# A piece is said by what it was carried on with, in the units pages use.

def test_a_piece_is_said_without_what_carrying_on_writes_into_every_piece(study):
    from fastmdxplora.gui.study_timeline import timeline_of

    piece = study / "segment-001"
    config = (study / "resolved_config.yml").read_text(encoding="utf-8")
    (piece / "resolved_config.yml").write_text(config.replace(
        "duration_ns: 1.0",
        "duration_ns: 0.0672\n  resume_from: /elsewhere/simulation/checkpoint.chk\n"
        "  setup_from: /elsewhere\n  prepared_from: /elsewhere\n  resume_unsealed: true\n"
        "  minimize: false\n  nvt_steps: 0\n  npt_steps: 0"), encoding="utf-8")
    said = next(e for e in timeline_of(study)["events"] if e["kind"] == "piece")
    assert [d["setting"] for d in said["settings"]] == ["simulation.duration_ns"]
    assert (said["settings"][0]["first_said"], said["settings"][0]["second_said"]) == \
        ("1 ns", "67.2 ps")


def test_the_side_panel_reads_again_after_a_switch_and_says_times_as_the_overview(study):
    """A switch while the Timeline's read was out drew the study before it,
    and nothing read the new one; times were on a 12-hour clock beside the
    Overview's 24, and lengths "30 min" beside its "30m 0s"."""
    pytest.importorskip("playwright.sync_api")
    from playwright.sync_api import sync_playwright

    from fastmdxplora.gui.server import start_dashboard_session

    session = start_dashboard_session(output=str(study), host="127.0.0.1", port=0)
    try:
        with sync_playwright() as pw:
            browser = pw.chromium.launch(args=["--enable-unsafe-swiftshader"])
            page = browser.new_page(viewport={"width": 1600, "height": 900}, locale="en-US")
            page.set_default_timeout(60000)
            errors: list[str] = []
            page.on("pageerror", lambda error: errors.append(str(error)))
            page.goto(session.url + "#overview", wait_until="domcontentloaded")
            page.wait_for_selector("body.state-ready")
            page.evaluate("""() => {
                const real = window.fetch;
                window.__asked = 0;
                window.fetch = (url, options) => {
                    if (!String(url).includes('/api/study-timeline')) return real(url, options);
                    window.__asked += 1;
                    if (window.__asked > 1) return real(url, options);
                    // The study before the switch, as its answer would say it.
                    return new Promise((resolve) => {
                        window.__release = () => resolve(new Response(JSON.stringify({
                            ok: true, events: [{kind: 'phase', phase: 'setup', said: 'STALE',
                                                at: '2026-10-18T07:00:00+00:00'}]})));
                    });
                };
                window.__drewStale = false;
                new MutationObserver(() => {
                    const list = document.getElementById('study-timeline-list');
                    if (list && list.textContent.includes('STALE')) window.__drewStale = true;
                }).observe(document.body, {childList: true, subtree: true});
            }""")
            page.evaluate("() => window.FastMDXFrame.setCollapsed(false)")
            page.click("#timeline-tab")
            page.wait_for_function("() => typeof window.__release === 'function'")
            # The study switches while the read is out.
            page.evaluate("""() => {
                const app = window.FastMDXDashboard.state.appState;
                window.FastMDXDashboard.applyAppState({...app, active_run: '/another/study'});
                window.__release();
            }""")
            page.wait_for_function("() => window.__asked === 2")
            page.wait_for_selector("#study-timeline-list .tl-event")
            page.click('.tl-event[data-kind="piece"] .tl-changes > summary')
            said = page.evaluate("""() => ({
                asked: window.__asked, stale: window.__drewStale,
                when: [...document.querySelectorAll('.tl-when')].map((n) => n.textContent),
                facts: document.querySelector('.tl-event[data-phase="simulation"] .tl-facts')
                    .textContent,
                heads: [...document.querySelectorAll('.tl-event[data-kind="piece"] .tl-table th')]
                    .map((n) => n.textContent),
            })""")
            browser.close()
    finally:
        session.server.shutdown()
    assert errors == []
    # The answer for the study before the switch is not drawn; the one open
    # is read again.
    assert said["stale"] is False and said["asked"] >= 2
    assert not any("AM" in w or "PM" in w for w in said["when"]), said["when"]
    assert said["facts"].startswith("20m 0s"), said["facts"]
    assert said["heads"] == ["Setting", "The study", "This piece"]


# ---------------------------------------------------------------------------
# A study of several runs takes no steps of its own: its runs' steps are
# listed, each named by its run, never "nothing recorded" when they ran.

def test_a_study_of_runs_lists_its_runs_steps_each_named_by_its_run(tmp_path):
    pytest.importorskip("playwright.sync_api")
    from playwright.sync_api import sync_playwright

    from fastmdxplora.gui.server import start_dashboard_session
    from fastmdxplora.gui.study_timeline import timeline_of
    from tests.test_the_overview_of_a_study_of_runs import _a_study_of_runs

    root = _a_study_of_runs(tmp_path, "simulation.random_seed", [1, 2, 3])
    runs = sorted((root / "runs").iterdir())
    for i, run in enumerate(runs[:2]):
        (run / "manifest.json").write_text(json.dumps({"phases": [
            _phase("setup", f"0{7 + i}:00:00", f"0{7 + i}:00:04"),
            _phase("simulation", f"0{7 + i}:00:04", f"0{7 + i}:30:04")]}), encoding="utf-8")
    (runs[2] / "manifest.json").unlink()

    said = timeline_of(root)
    assert said["ok"] is True
    assert [(e["run"], e["said"]) for e in said["events"]] == [
        ("seed 2", "Simulated"), ("seed 2", "Set up"), ("seed 1", "Simulated"), ("seed 1", "Set up")]
    assert [(r["label"], r["steps"]) for r in said["runs"]] == [
        ("seed 1", 2), ("seed 2", 2), ("seed 3", 0)]

    session = start_dashboard_session(output=str(root), host="127.0.0.1", port=0)
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
            shown = page.evaluate("""() => ({
                said: document.getElementById('study-timeline-said').textContent,
                steps: [...document.querySelectorAll('.tl-event .tl-said')].map((n) => n.textContent),
            })""")
            browser.close()
    finally:
        session.server.shutdown()
    assert errors == []
    assert "Nothing recorded" not in shown["said"]
    assert shown["said"] == ("The steps of its 3 runs, each named by the run that took it; "
                             "choose a run to see its own. Not yet recorded a step: seed 3.")
    assert shown["steps"] == ["seed 2: Simulated", "seed 2: Set up",
                              "seed 1: Simulated", "seed 1: Set up"]


# ---------------------------------------------------------------------------
# A step whose record gives no time says none, and the steps stay newest
# first by what is known.

def _timeless(run: Path, *names: str) -> None:
    """A run's Manifest whose phases record their name and status only, as
    a study of runs' runs may (no started_at, no finished_at)."""
    (run / "manifest.json").write_text(json.dumps({"phases": [
        {"name": name, "status": "ok"} for name in names]}), encoding="utf-8")


def test_a_step_with_no_time_says_none_and_comes_after_those_with_one(tmp_path):
    """Every step of a run whose Manifest gave no times read "Dec 31, 16:00"
    (the epoch on this computer's clock), and they came oldest first."""
    pytest.importorskip("playwright.sync_api")
    from playwright.sync_api import sync_playwright

    from fastmdxplora.gui.server import start_dashboard_session
    from fastmdxplora.gui.study_timeline import timeline_of
    from tests.test_the_overview_of_a_study_of_runs import _a_study_of_runs

    root = _a_study_of_runs(tmp_path, "simulation.random_seed", [1, 2])
    first, second = sorted((root / "runs").iterdir())
    _timeless(first, "setup", "simulation", "analysis", "report")
    (second / "manifest.json").write_text(json.dumps({"phases": [
        _phase("setup", "07:00:00", "07:00:04"), {"name": "simulation", "status": "ok"}]}),
        encoding="utf-8")

    said = timeline_of(root)
    # The one step with a time first; then those without, run by run in the
    # study's order, each run's newest first by its Manifest's order.
    assert [(e["run"], e["said"]) for e in said["events"]] == [
        ("seed 2", "Set up"),
        ("seed 1", "Report written"), ("seed 1", "Analysed"), ("seed 1", "Simulated"),
        ("seed 1", "Set up"), ("seed 2", "Simulated")]

    session = start_dashboard_session(output=str(root), host="127.0.0.1", port=0)
    try:
        with sync_playwright() as pw:
            browser = pw.chromium.launch(args=["--enable-unsafe-swiftshader"])
            page = browser.new_page(viewport={"width": 1600, "height": 900}, locale="en-US",
                                    timezone_id="America/Los_Angeles")
            page.set_default_timeout(60000)
            errors: list[str] = []
            page.on("pageerror", lambda error: errors.append(str(error)))
            page.goto(session.url + "#overview", wait_until="domcontentloaded")
            page.wait_for_selector("body.state-ready")
            page.evaluate("() => window.FastMDXFrame.setCollapsed(false)")
            page.click("#timeline-tab")
            page.wait_for_selector("#study-timeline-list .tl-event")
            shown = page.evaluate("""() => [...document.querySelectorAll('.tl-event')].map((n) => [
                n.querySelector('.tl-said').textContent,
                (n.querySelector('.tl-when') || {textContent: ''}).textContent])""")
            browser.close()
    finally:
        session.server.shutdown()
    assert errors == []
    assert [step for step, _ in shown] == [
        "seed 2: Set up", "seed 1: Report written", "seed 1: Analysed", "seed 1: Simulated",
        "seed 1: Set up", "seed 2: Simulated"]
    # The one time recorded is said; no step is given a time it has not.
    assert shown[0][1] != "" and all(when == "" for _, when in shown[1:]), shown


def test_a_run_asked_to_stop_is_said_stopped_as_its_overview_says_it(tmp_path):
    """A run stopped by Ctrl+C records its simulation as "error" with the
    refusal `simulation.run.stopped`: the Timeline said "Simulation error"
    where the Overview, the runs table and Compare said stopped."""
    from fastmdxplora.gui.overview_view import _phases
    from fastmdxplora.gui.study_timeline import timeline_of

    root = tmp_path / "study"
    root.mkdir()
    stopped = {**_phase("simulation", "07:00:04", "07:01:19"), "status": "error",
               "refusal": {"code": "simulation.run.stopped",
                           "message": "Production was stopped by Ctrl+C at step 6,400."}}
    (root / "manifest.json").write_text(json.dumps({"phases": [
        _phase("setup", "07:00:00", "07:00:04"), stopped]}), encoding="utf-8")
    said = {e["phase"]: e for e in timeline_of(root)["events"]}
    assert said["simulation"]["said"] == "Simulation stopped"
    assert said["simulation"]["status"] == {p["name"]: p for p in _phases(root)}[
        "simulation"]["status"] == "stopped"


# ---------------------------------------------------------------------------
# The rounds of a study run until what it asks is determined: each said by
# what was decided on it and, stopped, by what stopped it; placed by the
# recorded times of the steps around it, never below the oldest step.

def _round(decision: str, production: float, more: float | None = None,
           verdicts: list | None = None) -> dict:
    entry = {"production_ns": production, "decision": decision,
             "verdicts": verdicts or [{"analysis": "rmsd", "met": False,
                                       "said": "rmsd 0.150 ± 0.004 nm, outside the ±0.002 asked"}]}
    if more is not None:
        entry["more_ns"] = more
    return entry


def _stopping(where: Path, runs: list[str], rounds: list[dict], outcome: str, said: str) -> None:
    (where / "stopping.json").write_text(json.dumps({
        "targets": [{"analysis": "rmsd", "standard_error": 0.002}], "max_duration_ns": 5.0,
        "independent_starts": "required", "runs": runs, "rounds": rounds,
        "outcome": outcome, "said": said}), encoding="utf-8")


def _a_piece(run: Path, number: int, length_ns: float, start: str | None,
             end: str | None) -> None:
    piece = run / f"segment-{number:03d}"
    (piece / "simulation").mkdir(parents=True)
    config = (run / "resolved_config.yml").read_text(encoding="utf-8")
    (piece / "resolved_config.yml").write_text(
        config.replace("duration_ns: 1.0", f"duration_ns: {length_ns}"), encoding="utf-8")
    if start is not None:
        (piece / "manifest.json").write_text(json.dumps({"phases": [
            _phase("simulation", start, end)]}), encoding="utf-8")


@pytest.fixture()
def in_rounds(tmp_path) -> Path:
    """Two replicas run until RMSD is determined: round 1 extended each by
    500 ps, round 2 asked for 250 ps more and was stopped by Ctrl+C while
    the second run's piece ran."""
    from tests.test_the_overview_of_a_study_of_runs import _a_study_of_runs

    root = _a_study_of_runs(tmp_path, "simulation.random_seed", [1, 2])
    first, second = sorted((root / "runs").iterdir())
    for run, ready, end in ((first, "06:59:40", "07:30:00"), (second, "06:59:50", "07:40:00")):
        (run / "manifest.json").write_text(json.dumps({"phases": [
            _phase("setup", "06:59:00", ready),
            _phase("simulation", "07:00:00", end)]}), encoding="utf-8")
    _a_piece(first, 1, 0.5, "08:00:00", "08:20:00")
    _a_piece(second, 1, 0.5, "08:20:00", "08:40:00")
    _a_piece(first, 2, 0.25, "09:00:00", "09:10:00")
    _a_piece(second, 2, 0.25, None, None)
    _stopping(root, [first.name, second.name],
              [_round("extend", 1.0, 0.5), _round("stopped", 1.5, 0.25)], "stopped",
              f"Stopped when asked, while extending {second.name}. Resuming the study "
              "judges the runs again and carries on from there.")
    return root


def test_the_rounds_are_placed_by_when_they_happened_and_said_by_what_stopped_them(in_rounds):
    """The rounds came below the oldest step, and a round stopped by Ctrl+C
    was said "an extension failed"."""
    from fastmdxplora.gui.study_timeline import timeline_of

    said = timeline_of(in_rounds)
    steps = [(e.get("run", ""), e["said"]) for e in said["events"]]
    assert steps == [
        ("seed 1", "Carried on, piece 2"),
        ("", "Round 2: stopped when asked (Ctrl+C or a signal) while the runs were extended"),
        ("seed 2", "Carried on, piece 1"), ("seed 1", "Carried on, piece 1"),
        ("", "Round 1: each run extended by 500 ps"),
        ("seed 2", "Simulated"), ("seed 1", "Simulated"),
        ("seed 2", "Set up"), ("seed 1", "Set up")], steps
    rounds = [e for e in said["events"] if e["kind"] == "round"]
    assert [r["placed"] for r in rounds] == [
        {"at": "2026-10-18T09:00:00+00:00", "before": True},
        {"at": "2026-10-18T08:00:00+00:00", "before": True}]
    assert [r["production_said"] for r in rounds] == ["1.5 ns", "1 ns"]
    assert rounds[0]["stopped_by"] == "asked" and rounds[0]["at"] is None
    assert rounds[0]["detail"].startswith("Stopped when asked, while extending")
    assert not any("extension failed" in e["said"] for e in said["events"])


ROUND_STOPPED = pytest.mark.parametrize("rounds, outcome, said, last, by", [
    ([_round("stopped", 1.0, verdicts=[{"analysis": "rmsd", "met": False, "unrecorded": True,
                                        "said": "rmsd recorded no mean"}])],
     "stopped", "Stopped: rmsd recorded no mean. More production does not give an analysis "
     "a mean it does not record; its own log says why.",
     "Round 1: stopped: rmsd recorded no mean", "no mean"),
    ([_round("stopped", 1.0, 0.5)], "stopped",
     "Stopped: extending study failed at simulation: CUDA_ERROR_OUT_OF_MEMORY",
     "Round 1: stopped by an error while the runs were extended", "error"),
    # Stopped, then resumed: the record's sentence is the last round's, so
    # what stopped the first is not in it.
    ([_round("stopped", 1.0, 0.5), _round("met", 1.0)], "met",
     "After 1 ns of production: rmsd 0.150 ± 0.001 nm, within the ±0.002 asked.",
     "Round 1: stopped; what stopped it is not recorded", "not recorded"),
    ([_round("met", 1.0)], "met", "After 1 ns of production: rmsd 0.150 ± 0.001 nm.",
     "Round 1: determined as asked, so stopped", ""),
    ([_round("ceiling", 5.0)], "ceiling", "Stopped at the 5 ns ceiling.",
     "Round 1: stopped at the ceiling, not determined as asked", ""),
], ids=["no mean", "error", "not recorded", "met", "ceiling"])


@ROUND_STOPPED
def test_a_round_is_said_by_what_the_record_says_stopped_it(tmp_path, rounds, outcome, said,
                                                            last, by):
    from fastmdxplora.gui.study_timeline import timeline_of

    root = tmp_path / "study"
    root.mkdir()
    (root / "manifest.json").write_text(json.dumps({"phases": [
        _phase("setup", "07:00:00", "07:00:04"), _phase("simulation", "07:00:04", "07:30:04")]}),
        encoding="utf-8")
    _stopping(root, ["study"], rounds, outcome, said)
    events = timeline_of(root)["events"]
    first = next(e for e in events if e["kind"] == "round" and e["round"] == 1)
    assert (first["said"], first["stopped_by"]) == (last, by)
    # Placed after the production it judged, its own time not recorded.
    assert first["placed"] == {"at": "2026-10-18T07:30:04+00:00", "before": False}
    assert events[-2:] == [e for e in events if e["phase"] in ("simulation", "setup")]


def test_the_side_panel_says_a_round_s_length_in_ps_and_its_time_not_recorded(in_rounds):
    pytest.importorskip("playwright.sync_api")
    from playwright.sync_api import sync_playwright

    from fastmdxplora.gui.server import start_dashboard_session

    session = start_dashboard_session(output=str(in_rounds), host="127.0.0.1", port=0)
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
            page.wait_for_selector('#study-timeline-list .tl-event[data-kind="round"]')
            shown = page.evaluate("""() => [...document.querySelectorAll(
                '.tl-event[data-kind="round"]')].map((n) => ({
                said: n.querySelector('.tl-said').textContent,
                when: !!n.querySelector('.tl-when'),
                facts: n.querySelector('.tl-facts').textContent,
                by: n.dataset.stoppedBy || ''}))""")
            browser.close()
    finally:
        session.server.shutdown()
    assert errors == []
    assert [s["said"] for s in shown] == [
        "Round 2: stopped when asked (Ctrl+C or a signal) while the runs were extended",
        "Round 1: each run extended by 500 ps"]
    assert [s["by"] for s in shown] == ["asked", ""]
    assert not any(s["when"] for s in shown)
    assert shown[1]["facts"].startswith(
        "judged after 1 ns of production · its time not recorded: placed just before the "
        "pieces it began"), shown[1]["facts"]
