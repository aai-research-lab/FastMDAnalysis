"""A run the Agent started opens its own study; any other study opened
meanwhile is said as opened elsewhere (CI-774, fourth round).

Found by reviews K5 and K6: the first change of study in the minute after a
launch was taken as the run's own, whichever study it was, and a refused
launch left that minute open; a page that never saw the run's study kept
the run's thread and wrote in it; and with two launches, the first run's
study judged emptied the second's.
"""

from __future__ import annotations

import json
import time
import urllib.request

import pytest

from tests.test_a_conversation_follows_a_run_that_wrote_nothing import (  # noqa: F401
    _CONFIG,
    _a_study_opened_again_as_its_run_starts,
    _browsers_closed,
    _gui,
    _launches,
    _page_on,
    _say,
)

#: Each thing the page's notice says, in order, as it is shown.
_HEARD = """() => {
    window.__said = [];
    const note = () => {
        const toast = document.getElementById('dashboard-toast');
        const said = toast && toast.querySelector('.dashboard-toast-said');
        const text = toast && toast.classList.contains('show') && said ? said.textContent : '';
        if (text && window.__said[window.__said.length - 1] !== text) window.__said.push(text);
    };
    new MutationObserver(note).observe(document.body, {subtree: true, childList: true,
        characterData: true, attributes: true});
}"""

_RUN = "#agent-thread .agent-study:not([hidden]) [data-role=run]"


def _switch_elsewhere(session, folder):
    """Another window (or an AI app) opens ``folder``."""
    base = session.url.split("#")[0].rstrip("/")
    asked = urllib.request.Request(base + "/api/explore/switch", method="POST",
                                   data=json.dumps({"folder": str(folder)}).encode(),
                                   headers={"content-type": "application/json"})
    with urllib.request.urlopen(asked, timeout=30) as answer:
        assert json.loads(answer.read()).get("ok") is not False


def _study(folder):
    (folder / "simulation").mkdir(parents=True)
    (folder / "exploration.yml").write_text("systems: []\n", encoding="utf-8")
    return folder


def _shows(page, key, timeout=30000):
    page.wait_for_function("(key) => window.FastMDXDashboard.state.studyShown === key",
                           arg=key, timeout=timeout)


def _unsent(page):
    return page.evaluate("""() => {
        const kept = {};
        for (let i = 0; i < sessionStorage.length; i++) {
            const key = sessionStorage.key(i);
            if (key.startsWith('fmx.agentUnsent:')) {
                kept[key.slice('fmx.agentUnsent:'.length)] = JSON.parse(sessionStorage.getItem(key)).words;
            }
        }
        return kept;
    }""")


def _elsewhere_said(page):
    return [s for s in page.evaluate("() => window.__said") if "opened elsewhere" in s]


def test_a_run_s_study_never_seen_is_said_and_not_taken_as_this_window_s(tmp_path, monkeypatch):
    """K6 M1: another window opens the study the run was started from again
    as the run starts. The page says where the run's thread is; the run's
    study, opened later, is said as opened elsewhere and takes no words
    typed here."""
    from fastmdxplora.gui.workspace import study_key

    more = "typed after the run's thread went"

    def then(page, session, workspace):
        here, run = workspace / "here", workspace / "run1"
        page.evaluate(_HEARD)
        page.fill("#agent-request", more)
        page.wait_for_timeout(500)
        _switch_elsewhere(session, run)
        _shows(page, study_key(run))
        page.wait_for_function("() => window.__said.some((s) => /opened elsewhere/.test(s))",
                               timeout=15000)
        page.wait_for_timeout(1000)
        return {"said": _elsewhere_said(page), "box": page.input_value("#agent-request"),
                "kept": _unsent(page), "here": study_key(here), "run": study_key(run)}

    shown, workspace, words = _a_study_opened_again_as_its_run_starts(tmp_path, monkeypatch,
                                                                       then=then)
    assert "The run's conversation is in the run's study" in shown["note"], shown["note"]
    after = shown["then"]
    assert len(after["said"]) == 1, after["said"]
    assert more not in after["box"]
    assert after["kept"].get(after["here"]) == more, after["kept"]
    assert more not in after["kept"].get(after["run"], ""), after["kept"]
    assert shown["errors"] == []


@pytest.mark.parametrize("how", ["seen launching", "refused", "own"])
def test_a_study_opened_elsewhere_after_a_launch_is_said(tmp_path, monkeypatch, how):
    """K6 L1 and K5 2: a study another window opens while a run starts, and
    seen before the run's, is said as opened elsewhere; so is one opened
    after a launch that was refused. The run's own study is not."""
    from fastmdxplora.gui.server import start_dashboard_session
    from fastmdxplora.gui.workspace import study_key

    pytest.importorskip("playwright.sync_api")
    import tempfile

    monkeypatch.setenv("FASTMDXPLORA_CONFIG_DIR", tempfile.mkdtemp())
    workspace = tmp_path / "workspace"
    here = _study(workspace / "here")
    other = _study(workspace / "other")
    made = _launches(monkeypatch, workspace, taking=4.0)
    session = start_dashboard_session(output=str(here), host="127.0.0.1", port=0)
    held, answered = [], []
    page, errors, close = _page_on(session, [dict(_CONFIG)])
    try:
        page.on("response", lambda r: answered.append(r) if "/api/agent/run" in r.url else None)
        if how == "refused":
            page.route("**/api/agent/run", lambda route: route.fulfill(
                status=200, content_type="application/json",
                body=json.dumps({"ok": False, "error": "Not started here."})))
        _say(page, "chignolin")
        page.wait_for_selector(_RUN)
        _shows(page, study_key(here))
        page.wait_for_timeout(800)
        page.evaluate(_HEARD)
        if how == "seen launching":
            # The page reads no status while the GUI is on the run's folder.
            page.route("**/api/status*", lambda route: held.append(route))
        page.click(_RUN)
        if how == "seen launching":
            time.sleep(1.0)
            _switch_elsewhere(session, other)
            page.unroute("**/api/status*")
            for route in list(held):
                route.continue_()
            _shows(page, study_key(other))
        while not answered:
            page.wait_for_timeout(100)
        if how == "refused":
            page.wait_for_timeout(500)
            _switch_elsewhere(session, other)
            _shows(page, study_key(other))
        if how == "own":
            _shows(page, study_key(workspace / "run1"))
            page.wait_for_timeout(3000)
        else:
            page.wait_for_function(
                "() => window.__said.some((s) => /opened elsewhere/.test(s))", timeout=15000)
            page.wait_for_timeout(1000)
        said = _elsewhere_said(page)
        title = page.text_content("#topbar-run-title")
    finally:
        close()
        session.server.shutdown()
    if how == "own":
        assert said == [] and made
    else:
        assert said == [f"{title} was opened elsewhere, in another window or by an AI app: "
                        "this window now shows it."], said
    assert errors == []


def test_words_follow_the_thread_to_the_second_of_two_runs(tmp_path, monkeypatch):
    """K5 1: run A answered and not yet seen, run B started; the page sees
    A while B starts and B once it has answered. The words typed went with
    the thread to B, and are kept for B."""
    from fastmdxplora.gui.workspace import study_key

    workspace, made, session = _gui(tmp_path, monkeypatch, taking=0.5)
    words = "words for the thread"
    statuses, runs, later = [], [], []
    page, errors, close = _page_on(session, [dict(_CONFIG), dict(_CONFIG), dict(_CONFIG)])
    try:
        _say(page, "chignolin")
        page.wait_for_selector(_RUN)
        page.wait_for_timeout(800)
        page.click(_RUN)
        run1, run2, run3 = (study_key(workspace / f"run{i}") for i in (1, 2, 3))
        _shows(page, run1)
        page.wait_for_selector("#agent-thread .agent-running")
        _say(page, "again")
        page.wait_for_function(
            "() => document.querySelectorAll('#agent-thread .agent-study:not([hidden])').length === 2")
        page.wait_for_timeout(800)
        # The page does not see run A (run2) until B (run3) is launching.
        page.route("**/api/status*", lambda route: statuses.append(route))
        page.click(f"{_RUN}:visible >> nth=-1")
        page.wait_for_function(
            "() => document.querySelectorAll('#agent-thread .agent-running').length >= 2",
            timeout=30000)
        _say(page, "and again")
        page.wait_for_function(
            "() => document.querySelectorAll('#agent-thread .agent-study:not([hidden])').length === 3")
        page.wait_for_timeout(800)
        assert page.evaluate("() => window.FastMDXDashboard.state.studyShown") == run1
        page.fill("#agent-request", words)
        page.wait_for_timeout(500)
        page.route("**/api/agent/run", lambda route: runs.append(route))
        page.click(f"{_RUN}:visible >> nth=-1")
        while not runs:
            page.wait_for_timeout(100)
        page.unroute("**/api/status*")
        for route in list(statuses):
            route.continue_()
        _shows(page, run2, 15000)
        page.route("**/api/status*", lambda route: later.append(route))
        page.unroute("**/api/agent/run")
        for route in list(runs):
            route.continue_()
        page.wait_for_function(
            "() => document.querySelectorAll('#agent-thread .agent-running').length >= 3",
            timeout=30000)
        page.wait_for_timeout(300)
        page.unroute("**/api/status*")
        for route in list(later):
            route.continue_()
        _shows(page, run3)
        page.wait_for_timeout(1500)
        box = page.input_value("#agent-request")
        kept = _unsent(page)
        note = page.evaluate(
            "() => (document.getElementById('agent-unsent-note') || {}).textContent || ''")
    finally:
        close()
        session.server.shutdown()
    assert len(made) == 3
    assert words in box and kept.get(run3) == words, (box, kept, note)
    assert errors == []
