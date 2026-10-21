"""Notices saved, behind a bell in the status bar (IMDE 2.8).

A notice was a toast for three seconds and then gone; a run that ended
while the person was elsewhere said so only in the tab's title, and a
refusal only in the Log. Each is saved now on this computer, the last 200:
every toast as it is shown, whoever showed it; each run seen ending; each
refusal written while the page is open. The bell in the status bar says how
many are new; its list has All, Runs and Refusals and warnings, each
notice its time and its study, a run's end opening that study's Overview;
Clear the list asks nothing and offers Undo.
"""

from __future__ import annotations

import json
from datetime import datetime, timedelta, timezone

import pytest

LOOK = """() => {
    const bell = document.getElementById('sb-notices');
    const box = document.getElementById('notices');
    const count = bell.querySelector('.sb-count');
    return {
        bell: !bell.hidden, count: count.hidden ? '' : count.textContent,
        label: bell.getAttribute('aria-label'),
        open: !!box && !box.hidden,
        rows: box ? [...box.querySelectorAll('.notice')].map((n) => [n.dataset.kind,
            n.querySelector('.notice-said').textContent,
            (n.querySelector('.notice-act') || {}).textContent || '']) : [],
        foot: ((box && box.querySelector('.notices-foot')) || {}).textContent || '',
        focus: document.activeElement.id || document.activeElement.className,
        saved: JSON.parse(localStorage.getItem('fmx.notices') || '[]').length,
    };
}"""


def test_toasts_runs_ended_and_refusals_are_saved_behind_the_bell(tmp_path):
    pytest.importorskip("playwright.sync_api")
    from playwright.sync_api import sync_playwright

    from fastmdxplora.gui.server import start_dashboard_session
    from tests.test_the_workspace_says_its_studies import _study

    study = _study(tmp_path / "chignolin")
    session = start_dashboard_session(output=str(study), host="127.0.0.1", port=0)
    later = (datetime.now(timezone.utc) + timedelta(minutes=5)).isoformat()
    try:
        with sync_playwright() as pw:
            browser = pw.chromium.launch(args=["--enable-unsafe-swiftshader"])
            page = browser.new_page(viewport={"width": 1440, "height": 900})
            page.set_default_timeout(60000)
            errors: list[str] = []
            page.on("pageerror", lambda error: errors.append(str(error)))
            page.goto(session.url + "#overview", wait_until="domcontentloaded")
            page.wait_for_selector("body.state-ready")
            page.wait_for_function("() => !!window.FastMDXNotices")
            none = page.evaluate(LOOK)
            # A toast, from whichever page shows it; the same one again at
            # once is one notice.
            page.evaluate("() => window.FastMDXDashboard.toast('Saved the view.', 'ok')")
            page.evaluate("() => window.FastMDXDashboard.toast('The study could not open.', 'warning')")
            page.evaluate("() => window.FastMDXDashboard.toast('The study could not open.', 'warning')")
            # A run seen ending (run-notice.js), and a refusal written since
            # the Log was first read (frame.js).
            page.evaluate("""(study) => window.dispatchEvent(new CustomEvent('fmx:run-ended', {
                detail: {name: 'chignolin seed 2', run: 'completed', stage: 'production done',
                         study: study}}))""", str(study))
            page.evaluate("""() => window.dispatchEvent(new CustomEvent('fmx:events', {
                detail: {events: [], refused: () => true}}))""")
            page.evaluate("""(at) => window.dispatchEvent(new CustomEvent('fmx:events', {
                detail: {events: [{timestamp: at, message: 'simulation.timestep_fs 5 is refused'}],
                         refused: () => true}}))""", later)
            page.wait_for_function("() => window.FastMDXNotices.all.length === 4")
            gathered = page.evaluate(LOOK)
            page.click("#sb-notices")
            opened = page.evaluate(LOOK)
            page.click('.notices-chip[data-show="runs"]')
            runs = page.evaluate(LOOK)
            page.click('.notices-chip[data-show="refused"]')
            refused = page.evaluate(LOOK)
            page.click('.notices-chip[data-show="all"]')
            # Saved: there after a reload, read.
            page.reload(wait_until="domcontentloaded")
            page.wait_for_function("() => !!window.FastMDXNotices")
            reloaded = page.evaluate(LOOK)
            page.click("#sb-notices")
            page.click(".notices-clear")
            cleared = page.evaluate(LOOK)
            page.click(".notices-undo")
            undone = page.evaluate(LOOK)
            page.keyboard.press("Escape")
            closed = page.evaluate(LOOK)
            browser.close()
    finally:
        session.server.shutdown()
    assert errors == []
    assert none["bell"] and none["count"] == "" and none["label"] == "Notices"
    assert gathered["count"] == "4" and gathered["label"] == "Notices, 4 new" and gathered["saved"] == 4
    assert opened["open"] and opened["count"] == ""
    assert opened["rows"] == [
        ["refused", "simulation.timestep_fs 5 is refused", ""],
        ["ended", "chignolin seed 2 completed: production done", "Open its Overview"],
        ["warning", "The study could not open.", ""],
        ["ok", "Saved the view.", ""],
    ]
    assert [r[0] for r in runs["rows"]] == ["ended"]
    assert [r[0] for r in refused["rows"]] == ["refused", "warning"]
    assert reloaded["saved"] == 4 and reloaded["count"] == ""
    assert cleared["rows"] == [] and cleared["foot"].startswith("Cleared 4.") and cleared["saved"] == 0
    assert len(undone["rows"]) == 4 and undone["saved"] == 4
    assert undone["foot"].startswith("Saved on this computer, the last 200")
    assert not closed["open"] and closed["focus"] == "sb-notices"


def _session(tmp_path, *, allow_control=True):
    """Two studies of one workspace, the first open; served on loopback as
    a page bound beyond it where ``allow_control`` is False."""
    import contextlib
    import threading
    from http.server import ThreadingHTTPServer

    from fastmdxplora.gui.exploration import DashboardRuntime
    from fastmdxplora.gui.server import make_handler
    from tests.test_the_workspace_says_its_studies import _study

    root = (tmp_path / "workspace").resolve()
    alpha, beta = _study(root / "alpha"), _study(root / "beta")

    @contextlib.contextmanager
    def serving():
        runtime = DashboardRuntime(workspace_root=root, exploration_root=root, active_root=alpha)
        httpd = ThreadingHTTPServer(("127.0.0.1", 0),
                                    make_handler(alpha, runtime=runtime, allow_control=allow_control))
        threading.Thread(target=httpd.serve_forever, daemon=True).start()
        try:
            yield f"http://127.0.0.1:{httpd.server_address[1]}"
        finally:
            httpd.shutdown()
            httpd.server_close()
    return alpha, beta, serving()


ENDED = """(study) => window.dispatchEvent(new CustomEvent('fmx:run-ended', {
    detail: {name: 'beta', run: 'completed', stage: 'production done', study: study}}))"""


def test_a_notice_opens_another_study_only_as_every_switch_does(tmp_path):
    pytest.importorskip("playwright.sync_api")
    from playwright.sync_api import sync_playwright

    alpha, beta, serving = _session(tmp_path)
    with serving as url, sync_playwright() as pw:
        browser = pw.chromium.launch(args=["--enable-unsafe-swiftshader"])
        page = browser.new_page(viewport={"width": 1440, "height": 900})
        page.set_default_timeout(60000)
        errors: list[str] = []
        page.on("pageerror", lambda error: errors.append(str(error)))
        switches: list[str] = []
        page.on("request", lambda r: switches.append(r.url) if "/api/explore/switch" in r.url else None)
        page.goto(url + "#analysis", wait_until="domcontentloaded")
        page.wait_for_selector("body.state-ready")
        page.wait_for_function("() => !!window.FastMDXNotices")
        page.evaluate(ENDED, str(beta))
        # Work here a switch would lose: the notice's way to the other
        # study is held as the strip's and the sidebar's are.
        page.evaluate("""() => { window.FastMDXRun = window.FastMDXRun || {};
            window.FastMDXRun.beforeStudySwitch = () => 'The draft here is not saved.'; }""")
        page.click("#sb-notices")
        page.click(".notice-act")
        page.wait_for_function("""() => (document.querySelector('#dashboard-toast .dashboard-toast-said')
            || {}).textContent === 'The draft here is not saved.'""")
        held = {"switches": list(switches), "hash": page.evaluate("location.hash"),
                "open": page.evaluate("window.FastMDXDashboard.state.appState.active_run")}
        page.evaluate("() => { window.FastMDXRun.beforeStudySwitch = () => ''; }")
        page.click("#sb-notices")
        page.click(".notice-act")
        page.wait_for_function("(beta) => window.FastMDXDashboard.state.appState.active_run === beta",
                               arg=str(beta))
        page.wait_for_function("() => location.hash.startsWith('#overview')")
        browser.close()
    assert errors == []
    assert held == {"switches": [], "hash": "#analysis", "open": str(alpha)}
    assert len(switches) == 1


def test_a_notice_saves_only_what_it_said(tmp_path):
    pytest.importorskip("playwright.sync_api")
    from playwright.sync_api import sync_playwright

    alpha, beta, serving = _session(tmp_path)
    later = (datetime.now(timezone.utc) + timedelta(minutes=5)).isoformat()
    with serving as url, sync_playwright() as pw:
        browser = pw.chromium.launch(args=["--enable-unsafe-swiftshader"])
        page = browser.new_page(viewport={"width": 1440, "height": 900})
        page.set_default_timeout(60000)
        page.goto(url + "#overview", wait_until="domcontentloaded")
        page.wait_for_selector("body.state-ready")
        page.wait_for_function("() => !!window.FastMDXNotices")
        page.evaluate("() => window.FastMDXDashboard.toast('Saved the view.', 'ok')")
        page.evaluate("""() => window.dispatchEvent(new CustomEvent('fmx:events', {
            detail: {events: [], refused: () => true}}))""")
        page.evaluate("""([at, where]) => window.dispatchEvent(new CustomEvent('fmx:events', {
            detail: {events: [{timestamp: at, message: 'analysis.rmsd.selection is refused: no atoms\\n'
                + 'Traceback (most recent call last):\\n  File "' + where + '/secret_notes.py", line 3'}],
                     refused: () => true}}))""", [later, str(alpha)])
        page.evaluate(ENDED, str(beta))
        page.wait_for_function("() => window.FastMDXNotices.all.length === 3")
        saved = page.evaluate("() => JSON.parse(localStorage.getItem('fmx.notices'))")
        browser.close()
    said = {n["kind"]: n for n in saved}
    assert said["refused"]["said"] == "analysis.rmsd.selection is refused: no atoms"
    # The study's folder only where it opens the study (a run's end).
    assert said["ok"]["folder"] == "" and said["refused"]["folder"] == ""
    assert said["ended"]["folder"] == str(beta)
    assert "secret_notes" not in json.dumps(saved) and str(alpha) not in json.dumps(saved)
