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


def _status(study, said, step=5):
    """The run's live record, as a run writes it."""
    (study / "simulation").mkdir(exist_ok=True)
    (study / "simulation" / "live_status.json").write_text(json.dumps({
        "status": said, "stage": "production", "system_id": "TRIA", "current_step": step,
        "total_planned_steps": 10, "run_started_at": "2026-10-21T10:00:00+00:00",
        "last_update_timestamp": datetime.now(timezone.utc).isoformat()}), encoding="utf-8")


RUNNING = "() => document.getElementById('status-run').getAttribute('data-run') === 'running'"
SAVED = "() => JSON.parse(localStorage.getItem('fmx.notices') || '[]')"


def _page(browser, url, where="#analysis", width=1440, context=None):
    page = (context or browser.new_context(viewport={"width": width, "height": 900})).new_page()
    page.set_default_timeout(60000)
    page.errors = []
    page.on("pageerror", lambda error: page.errors.append(str(error)))
    page.goto(url + where, wait_until="domcontentloaded")
    page.wait_for_selector("body.state-ready")
    page.wait_for_function("() => !!window.FastMDXNotices")
    return page


def _run_ends(study, pages, said="completed"):
    """A run of ``study`` seen going in each page, then ended."""
    _status(study, "running")
    for page in pages:
        page.wait_for_function(RUNNING)
    _status(study, said, step=10 if said == "completed" else 5)
    for page in pages:
        page.wait_for_function(f"() => ({SAVED})().some((n) => n.act === 'overview')")


def test_a_run_s_end_opens_its_study_by_its_key(tmp_path):
    """A real run's end: saved with the server's key for its study and its
    folder's name; Open its Overview shows its own study's Overview, and
    another study's opened by its key as every switch is, held first."""
    pytest.importorskip("playwright.sync_api")
    from playwright.sync_api import sync_playwright

    from fastmdxplora.gui.workspace import study_key

    alpha, beta, serving = _session(tmp_path)
    with serving as url, sync_playwright() as pw:
        browser = pw.chromium.launch(args=["--enable-unsafe-swiftshader"])
        page = _page(browser, url)
        switches: list[str] = []
        page.on("request", lambda r: switches.append(r.post_data or "")
                if "/api/explore/switch" in r.url else None)
        _run_ends(alpha, [page])
        saved = page.evaluate(SAVED)
        page.click("#sb-notices")
        page.click(".notice-act")
        page.wait_for_function("() => location.hash.startsWith('#overview')")
        own = list(switches)
        # Another study open: the run's end opens alpha again, by its key.
        page.evaluate("(beta) => window.FastMDXPlace.open(beta, 'analysis')", str(beta))
        page.wait_for_function("(beta) => !!window.FastMDXDashboard && (window.FastMDXDashboard.state.appState || {}).active_run === beta",
                               arg=str(beta))
        page.wait_for_selector("body.state-ready")
        page.wait_for_function("() => !!window.FastMDXNotices")
        page.evaluate("""() => { window.FastMDXRun = window.FastMDXRun || {};
            window.FastMDXRun.beforeStudySwitch = () => 'The draft here is not saved.'; }""")
        before = len(switches)
        page.click("#sb-notices")
        page.click(".notice-act")
        page.wait_for_function("""() => (document.querySelector('#dashboard-toast .dashboard-toast-said')
            || {}).textContent === 'The draft here is not saved.'""")
        held = len(switches) - before
        page.evaluate("() => { window.FastMDXRun.beforeStudySwitch = () => ''; }")
        page.click("#sb-notices")
        page.click(".notice-act")
        page.wait_for_function("(alpha) => !!window.FastMDXDashboard && (window.FastMDXDashboard.state.appState || {}).active_run === alpha",
                               arg=str(alpha))
        page.wait_for_function("() => location.hash.startsWith('#overview')")
        errors = page.errors
        browser.close()
    assert errors == []
    ended = [n for n in saved if n["act"] == "overview"]
    assert len(ended) == 1 and ended[0]["kind"] == "ended"
    assert ended[0]["key"] == study_key(alpha) and ended[0]["name"] == "alpha"
    assert "folder" not in ended[0] and str(alpha) not in json.dumps(saved)
    assert own == [] and held == 0
    assert json.loads(switches[-1]) == {"study": study_key(alpha)}


def test_a_run_that_fails_with_no_record_is_saved_as_it_failed(tmp_path):
    """A run with no record of its own (a study of several runs, a process
    ended before it wrote one) ends as its process did: failed, under
    Refusals and warnings, as the browser's notice says it."""
    pytest.importorskip("playwright.sync_api")
    from playwright.sync_api import sync_playwright

    alpha, beta, serving = _session(tmp_path)
    process = {"state": "idle"}

    def app_state(route):
        answer = route.fetch().json()
        if process["state"] == "running":
            answer.update(process_running=True, status="running")
        elif process["state"] == "failed":
            answer.update(process_running=False, status="failed")
        route.fulfill(json=answer)

    with serving as url, sync_playwright() as pw:
        browser = pw.chromium.launch(args=["--enable-unsafe-swiftshader"])
        context = browser.new_context(viewport={"width": 1440, "height": 900})
        context.add_init_script("""window.__shown = [];
            window.Notification = function (title) { window.__shown.push(title); };
            window.Notification.permission = 'granted';""")
        context.route("**/api/app-state*", app_state)
        page = _page(browser, url, context=context)
        page.evaluate("() => { document.getElementById('setting-notify').checked = true; }")
        process["state"] = "running"
        page.wait_for_function(RUNNING)
        process["state"] = "failed"
        page.wait_for_function(f"() => ({SAVED})().some((n) => n.act === 'overview')")
        saved = [n for n in page.evaluate(SAVED) if n["act"] == "overview"]
        shown = page.evaluate("window.__shown")
        page.click("#sb-notices")
        page.click('.notices-chip[data-show="refused"]')
        refused = page.evaluate("[...document.querySelectorAll('.notice .notice-said')].map((n) => n.textContent)")
        browser.close()
    assert len(saved) == 1 and saved[0]["kind"] == "refused" and " failed" in saved[0]["said"]
    assert len(shown) == 1 and shown[0].endswith(" failed")
    assert saved[0]["said"] in refused


def test_a_notice_saves_only_what_it_said(tmp_path):
    pytest.importorskip("playwright.sync_api")
    from playwright.sync_api import sync_playwright

    alpha, beta, serving = _session(tmp_path)
    later = (datetime.now(timezone.utc) + timedelta(minutes=5)).isoformat()
    with serving as url, sync_playwright() as pw:
        browser = pw.chromium.launch(args=["--enable-unsafe-swiftshader"])
        page = _page(browser, url, "#overview")
        page.evaluate("() => window.FastMDXDashboard.toast('Saved the view.', 'ok')")
        page.evaluate("""() => window.dispatchEvent(new CustomEvent('fmx:events', {
            detail: {events: [], refused: () => true}}))""")
        page.evaluate("""([at, where]) => window.dispatchEvent(new CustomEvent('fmx:events', {
            detail: {events: [{timestamp: at, message: 'analysis.rmsd.selection is refused: no atoms\\n'
                + 'Traceback (most recent call last):\\n  File "' + where + '/secret_notes.py", line 3'}],
                     refused: () => true}}))""", [later, str(alpha)])
        page.wait_for_function("() => window.FastMDXNotices.all.length === 2")
        saved = page.evaluate(SAVED)
        browser.close()
    said = {n["kind"]: n for n in saved}
    assert said["refused"]["said"] == "analysis.rmsd.selection is refused: no atoms"
    # No study's folder: a key only where it opens the study (a run's end).
    assert said["ok"]["key"] == "" and said["refused"]["key"] == ""
    assert "secret_notes" not in json.dumps(saved) and str(alpha) not in json.dumps(saved)


def test_a_toast_is_filed_under_the_study_of_its_page(tmp_path):
    """A notice on a page of no study (All studies, Compare) is about
    another study or none: never filed under the open one."""
    pytest.importorskip("playwright.sync_api")
    from playwright.sync_api import sync_playwright

    alpha, beta, serving = _session(tmp_path)
    with serving as url, sync_playwright() as pw:
        browser = pw.chromium.launch(args=["--enable-unsafe-swiftshader"])
        page = _page(browser, url, "#overview")
        page.evaluate("() => window.FastMDXDashboard.toast('Saved the view.', 'ok')")
        page.evaluate("() => window.FastMDXDashboard.navigate('studies')")
        page.wait_for_function("() => document.documentElement.dataset.page === 'studies'")
        page.evaluate("() => window.FastMDXDashboard.toast('Saved the tags and note of beta.', 'ok')")
        page.wait_for_function("() => window.FastMDXNotices.all.length === 2")
        named = {n["said"]: n["name"] for n in page.evaluate(SAVED)}
        browser.close()
    assert named == {"Saved the view.": "alpha", "Saved the tags and note of beta.": ""}


READ_THE_LOG = """() => fetch('/api/app-state').then((r) => r.json()).then((d) =>
    window.dispatchEvent(new CustomEvent('dashboard:app-state', {detail: d})))"""


def test_two_windows_save_one_notice_once(tmp_path):
    """A refusal the Log reads and a run's end, seen by two windows of the
    GUI, are saved once."""
    pytest.importorskip("playwright.sync_api")
    from playwright.sync_api import sync_playwright

    alpha, beta, serving = _session(tmp_path)
    with serving as url, sync_playwright() as pw:
        browser = pw.chromium.launch(args=["--enable-unsafe-swiftshader"])
        context = browser.new_context(viewport={"width": 1440, "height": 900})
        first, second = _page(browser, url, context=context), _page(browser, url, context=context)
        _status(alpha, "running")
        for page in (first, second):
            page.wait_for_function(RUNNING)
        events = alpha / "simulation" / "live_events.log"
        events.write_text("", encoding="utf-8")
        for page in (first, second):
            # Each window has read the open study's Log once, as it does with
            # each app state: what is written after is new.
            with page.expect_response(lambda r: "/api/events" in r.url):
                page.evaluate(READ_THE_LOG)
        at = datetime.now(timezone.utc).isoformat()
        events.write_text(f"{at}\terror\tsetup.forcefield.unknown is refused\n"
                          f"{at}\terror\tFastMDXplora exploration failed\n", encoding="utf-8")
        for page in (first, second):
            # Each window reads the Log again, as it does with each app state.
            with page.expect_response(lambda r: "/api/events" in r.url):
                page.evaluate(READ_THE_LOG)
            page.wait_for_timeout(300)
        _run_ends(alpha, [], said="completed")
        for page in (first, second):
            page.wait_for_function(f"() => ({SAVED})().some((n) => n.act === 'overview')")
        first.wait_for_timeout(3000)
        saved = first.evaluate(SAVED)
        browser.close()
    kinds = sorted(n["kind"] for n in saved)
    assert kinds == ["ended", "refused", "refused"], saved


def test_keys_pressed_in_the_list_stay_in_it(tmp_path):
    pytest.importorskip("playwright.sync_api")
    from playwright.sync_api import sync_playwright

    alpha, beta, serving = _session(tmp_path)
    chosen = ("() => [...document.querySelectorAll('#studies-table tbody input[type=checkbox]')]"
              ".filter((b) => b.checked).length")
    with serving as url, sync_playwright() as pw:
        browser = pw.chromium.launch(args=["--enable-unsafe-swiftshader"])
        page = _page(browser, url, "#studies")
        page.evaluate(f"() => window.FastMDXStudies.load({str(alpha.parent)!r})")
        page.wait_for_selector("#studies-table tbody tr.studies-row")
        for box in page.query_selector_all("#studies-table tbody tr.studies-row input[type='checkbox']")[:2]:
            box.check()
        page.evaluate("() => window.FastMDXDashboard.toast('A notice.', 'ok')")
        page.click("#sb-notices")
        for key in ("c", "j", "k", "t", "x"):
            page.keyboard.press(key)
        typed = {"hash": page.evaluate("location.hash"), "chosen": page.evaluate(chosen),
                 "focus": page.evaluate("document.getElementById('notices').contains(document.activeElement)")}
        page.keyboard.press("Escape")
        closed = {"open": page.evaluate("!document.getElementById('notices').hidden"),
                  "chosen": page.evaluate(chosen), "focus": page.evaluate("document.activeElement.id")}
        browser.close()
    assert typed == {"hash": "#studies", "chosen": 2, "focus": True}
    assert closed == {"open": False, "chosen": 2, "focus": "sb-notices"}


def test_undo_after_clear_keeps_the_newest_first_and_the_count_right(tmp_path):
    pytest.importorskip("playwright.sync_api")
    from playwright.sync_api import sync_playwright

    alpha, beta, serving = _session(tmp_path)
    count = "() => { const c = document.querySelector('#sb-notices .sb-count'); return c.hidden ? '' : c.textContent; }"
    with serving as url, sync_playwright() as pw:
        browser = pw.chromium.launch(args=["--enable-unsafe-swiftshader"])
        page = _page(browser, url, "#overview")
        page.evaluate("() => window.FastMDXDashboard.toast('Older one.', 'ok')")
        page.wait_for_timeout(50)
        page.evaluate("() => window.FastMDXDashboard.toast('Older two.', 'ok')")
        page.click("#sb-notices")
        page.click(".notices-clear")
        page.wait_for_timeout(50)
        page.evaluate("() => window.FastMDXDashboard.toast('Newest, after Clear.', 'warning')")
        page.click(".notices-undo")
        order = [n["said"] for n in page.evaluate("window.FastMDXNotices.all")]
        page.keyboard.press("Escape")
        page.click("#sb-notices")
        page.keyboard.press("Escape")
        seen = page.evaluate(count)
        browser.close()
    assert order == ["Newest, after Clear.", "Older two.", "Older one."]
    assert seen == ""


def test_a_page_bound_beyond_this_computer_offers_no_switch_or_stop(tmp_path):
    pytest.importorskip("playwright.sync_api")
    from playwright.sync_api import sync_playwright

    from fastmdxplora.gui.workspace import study_key

    alpha, beta, serving = _session(tmp_path, allow_control=False)
    other = [{"at": 1, "kind": "ended", "said": "beta completed.", "name": "beta", "act": "overview",
              "key": study_key(beta)}]
    with serving as url, sync_playwright() as pw:
        browser = pw.chromium.launch(args=["--enable-unsafe-swiftshader"])
        context = browser.new_context(viewport={"width": 1440, "height": 900})
        # Another study's run's end, saved by another window.
        context.add_init_script(f"localStorage.setItem('fmx.notices', {json.dumps(json.dumps(other))});")
        page = _page(browser, url, "#overview", context=context)
        # Its own study's run, going, then ended.
        _status(alpha, "running")
        page.wait_for_function(RUNNING)
        stops = page.evaluate("[...document.querySelectorAll('[data-stop-run]')].filter((b) => !b.hidden).length")
        running = page.evaluate("window.FastMDXStop.running()")
        _status(alpha, "completed", step=10)
        page.wait_for_function(f"() => ({SAVED})().length === 2")
        page.click("#sb-notices")
        acts = page.evaluate("""[...document.querySelectorAll('.notice')].map((n) =>
            (n.querySelector('.notice-act') || {}).textContent || '')""")
        browser.close()
    # Its own study's Overview is a page of this one; another study is a
    # switch, a POST refused here.
    assert acts == ["Open its Overview", ""]
    assert stops == 0 and running is False
