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
import re
import time
import urllib.request

import pytest

from tests.test_a_conversation_follows_a_run_that_wrote_nothing import (  # noqa: F401
    _CONFIG,
    _a_study_opened_again_as_its_run_starts,
    _answer,
    _browsers_closed,
    _entries_everywhere,
    _gui,
    _launches,
    _page_on,
    _say,
)
from tests.test_open_studies_have_tabs import _builder, _open, browser, site  # noqa: F401

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
    """K5 1: run A answered and not yet seen, a message asking for run B is
    not sent (Send waits, saying so) until the page has read which study is
    open (it sees A) and the person sends it; B is started, and the words
    typed go with the thread to B, and are kept for B."""
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
        _when_known(page)
        _say(page, "again")
        page.wait_for_function(
            "() => document.querySelectorAll('#agent-thread .agent-study:not([hidden])').length === 2")
        page.wait_for_timeout(800)
        # The page does not see run A (run2) yet.
        page.route("**/api/status*", lambda route: statuses.append(route))
        page.click(f"{_RUN}:visible >> nth=-1")
        page.wait_for_function(
            "() => document.querySelectorAll('#agent-thread .agent-running').length >= 2",
            timeout=30000)
        _say(page, "and again")
        page.wait_for_timeout(1500)
        # Waits: the thread went to run A, which the page has not seen.
        held = _note(page)
        cards = page.locator("#agent-thread .agent-study:not([hidden])").count()
        assert page.evaluate("() => window.FastMDXDashboard.state.studyShown") == run1
        page.unroute("**/api/status*")
        for route in list(statuses):
            route.continue_()
        _shows(page, run2)
        _send_when_offered(page)
        page.wait_for_function(
            "() => document.querySelectorAll('#agent-thread .agent-study:not([hidden])').length === 3")
        page.wait_for_timeout(800)
        page.fill("#agent-request", words)
        page.wait_for_timeout(500)
        page.route("**/api/status*", lambda route: later.append(route))
        page.route("**/api/agent/run", lambda route: runs.append(route))
        page.click(f"{_RUN}:visible >> nth=-1")
        while not runs:
            page.wait_for_timeout(100)
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
        note = _note(page)
    finally:
        close()
        session.server.shutdown()
    assert len(made) == 3
    assert cards == 2 and held == _SEND_WAITS, (cards, held)
    assert words in box and kept.get(run3) == words, (box, kept, note)
    assert errors == []


def _here_and_other(tmp_path, monkeypatch, taking=4.0):
    import tempfile

    pytest.importorskip("playwright.sync_api")
    from fastmdxplora.gui.server import start_dashboard_session

    monkeypatch.setenv("FASTMDXPLORA_CONFIG_DIR", tempfile.mkdtemp())
    workspace = tmp_path / "workspace"
    here = _study(workspace / "here")
    other = _study(workspace / "other")
    made = _launches(monkeypatch, workspace, taking=taking)
    session = start_dashboard_session(output=str(here), host="127.0.0.1", port=0)
    return workspace, here, other, made, session


def _when_known(page):
    """Once the page has read which study is open after a launch: Send is
    offered then, not while the run starts (CI-774 round 3)."""
    page.wait_for_function("() => { const d = window.FastMDXDashboard;"
                           " return d.studyKnown() && d.viewKnown(); }", timeout=30000)
    page.wait_for_timeout(300)


def _note(page):
    return page.evaluate("() => (document.getElementById('agent-unsent-note') || {}).textContent || ''")


def test_a_study_seen_after_the_answer_leaves_no_run_s_study_this_window_s(tmp_path, monkeypatch):
    """L1, order A: another window opens `other` during the launch and the
    page sees it only after the answer; another window then opens the run's
    study. It is said as opened elsewhere, and the words typed for `other`
    stay kept for `other`."""
    from fastmdxplora.gui.workspace import study_key

    workspace, here, other, made, session = _here_and_other(tmp_path, monkeypatch)
    run = workspace / "run1"
    held, answered = [], []
    page, errors, close = _page_on(session, [dict(_CONFIG)])
    try:
        page.on("response", lambda r: answered.append(r) if "/api/agent/run" in r.url else None)
        _say(page, "chignolin")
        page.wait_for_selector(_RUN)
        _shows(page, study_key(here))
        page.wait_for_timeout(800)
        page.evaluate("""() => { window.__else = [];
            window.addEventListener('dashboard:study-changed-elsewhere',
                                    (e) => window.__else.push(e.detail.to)); }""")
        page.route("**/api/status*", lambda route: held.append(route))
        page.click(_RUN)
        page.focus("#agent-request")
        page.keyboard.type("typed during launch")
        time.sleep(1.0)
        _switch_elsewhere(session, other)
        while not answered:
            page.wait_for_timeout(100)
        page.wait_for_timeout(300)
        page.unroute("**/api/status*")
        for route in list(held):
            route.continue_()
        _shows(page, study_key(other))
        page.wait_for_timeout(2000)
        page.fill("#agent-request", "typed for other")
        page.wait_for_timeout(600)
        _switch_elsewhere(session, run)
        _shows(page, study_key(run))
        page.wait_for_function("() => window.__else.length >= 2", timeout=15000)
        page.wait_for_timeout(1000)
        said = page.evaluate("() => window.__else")
        box = page.input_value("#agent-request")
        kept = _unsent(page)
    finally:
        close()
        session.server.shutdown()
    assert made
    assert said == [study_key(other), study_key(run)], said
    assert "typed for other" not in box
    assert kept.get(study_key(other)) == "typed for other", kept
    assert "typed for other" not in kept.get(study_key(run), ""), kept
    assert errors == []


def test_the_first_of_two_runs_never_seen_is_said_when_opened_elsewhere(tmp_path, monkeypatch):
    """L1, order B: the Agent's run A (run1) answered and not yet seen, the
    Config Builder's run B (run2) started; the page sees only B. Another
    window then opens A: said as opened elsewhere, and the words typed in
    B are kept for B, not taken to A."""
    from fastmdxplora.gui.workspace import study_key

    workspace, here, other, made, session = _here_and_other(tmp_path, monkeypatch, taking=0.5)
    statuses, answered, built = [], [], []
    page, errors, close = _page_on(session, [dict(_CONFIG)])
    try:
        page.on("response", lambda r: answered.append(r) if "/api/agent/run" in r.url else None)
        page.on("response", lambda r: built.append(r) if r.url.endswith("/api/run") else None)
        _say(page, "chignolin")
        page.wait_for_selector(_RUN)
        _shows(page, study_key(here))
        page.wait_for_timeout(800)
        page.evaluate("""() => { window.__else = [];
            window.addEventListener('dashboard:study-changed-elsewhere',
                                    (e) => window.__else.push(e.detail.to)); }""")
        # The page does not see run A until B has started.
        page.route("**/api/status*", lambda route: statuses.append(route))
        page.click(_RUN)
        while not answered:
            page.wait_for_timeout(100)
        _builder(page, {"root": workspace})
        page.click("#run-start-button")
        while not built:
            page.wait_for_timeout(100)
        page.unroute("**/api/status*")
        for route in list(statuses):
            route.continue_()
        run_a, run_b = workspace / "run1", workspace / "run2"
        _shows(page, study_key(run_b))
        page.wait_for_timeout(1500)
        page.evaluate("window.FastMDXDashboard.navigate('agent')")
        page.fill("#agent-request", "typed in B")
        page.wait_for_timeout(600)
        _switch_elsewhere(session, run_a)
        _shows(page, study_key(run_a))
        page.wait_for_function("() => window.__else.length >= 1", timeout=15000)
        page.wait_for_timeout(1000)
        said = page.evaluate("() => window.__else")
        box = page.input_value("#agent-request")
        kept = _unsent(page)
    finally:
        close()
        session.server.shutdown()
    assert len(made) == 2
    assert said == [study_key(run_a)], said
    assert kept.get(study_key(run_b)) == "typed in B", kept
    assert "typed in B" not in box
    assert errors == []


@pytest.mark.parametrize("sent", ["during the launch", 0.3, 1.5])
def test_a_message_sent_as_the_run_starts_is_answered_in_the_study_shown(tmp_path, monkeypatch,
                                                                          sent):
    """L1 M2 and L2 M1: another window opens the start study again as the
    run starts, and the message is sent while the run starts or just after
    it answers, before the page's next status. Send waits, saying so; once
    the page knows, the person sends it, and it is answered with the study
    shown open and kept in that study's thread, never in the run's."""
    from fastmdxplora.gui.workspace import study_key

    workspace, here, other, made, session = _here_and_other(tmp_path, monkeypatch)
    words = "words typed as it starts"
    held, answered, asked = [], [], []
    page, errors, close = _page_on(session, [dict(_CONFIG), _answer("Said here.")])
    try:
        page.on("response", lambda r: answered.append(r) if "/api/agent/run" in r.url else None)
        page.on("request", lambda r: asked.append(json.loads(r.post_data or "{}"))
                if "/api/agent/propose" in r.url else None)
        _say(page, "chignolin")
        page.wait_for_selector(_RUN)
        _shows(page, study_key(here))
        page.wait_for_timeout(800)
        page.route("**/api/status*", lambda route: held.append(route))
        page.click(_RUN)
        page.focus("#agent-request")
        page.keyboard.type(words)
        time.sleep(1.0)
        _switch_elsewhere(session, here)
        page.unroute("**/api/status*")
        released = list(held)
        for route in released:
            route.continue_()
        if sent == "during the launch":
            page.wait_for_timeout(1200)
            assert not answered
        else:
            while not answered:
                page.wait_for_timeout(20)
            # No status read comes before the send.
            page.route("**/api/status*", lambda route: held.append(route))
            page.wait_for_timeout(int(sent * 1000))
        page.focus("#agent-request")
        page.keyboard.press("Enter")
        page.wait_for_timeout(300)
        waited = (_note(page) if page.input_value("#agent-request") else "sent at once") \
            if sent == "during the launch" else ""
        if sent != "during the launch":
            page.wait_for_timeout(1000)
            page.unroute("**/api/status*")
            for route in held[len(released):]:
                route.continue_()
        _send_when_offered(page)
        page.wait_for_selector("#agent-thread .agent-answer:has-text('Said here.')", timeout=30000)
        page.wait_for_timeout(1500)
        thread = page.inner_text("#agent-thread")
    finally:
        close()
        session.server.shutdown()
    kept = _entries_everywhere(workspace)
    in_here = [e.get("text") for path, entries in kept.items()
               if "here/agent" in path for e in entries]
    in_run = [e.get("text") for path, entries in kept.items()
              if "run1/agent" in path for e in entries]
    runs = [path.rsplit("/", 1)[1][:-5] for path in kept if "run1/agent" in path]
    assert made and runs
    assert words not in in_run and "Said here." not in in_run, kept
    assert words in in_here and "Said here." in in_here, kept
    assert "Said here." in thread and "chignolin" not in thread
    # While the launch is unanswered Send waits; after it, the message is
    # sent at once or when the person sends it again, into the study
    # shown's own thread only (asserted below).
    assert sent != "during the launch" or waited == _SEND_WAITS, waited
    assert [a.get("conversation") for a in asked if a.get("request") == words
            and a.get("conversation") in runs] == []
    assert errors == []


def test_words_typed_for_the_run_s_thread_follow_it_past_another_study(tmp_path, monkeypatch):
    """L2 L1: words typed after Run into the thread the run took with it,
    another study seen first: they are kept for the run's study, where the
    thread went, not for the study the run was started from."""
    from fastmdxplora.gui.workspace import study_key

    workspace, here, other, made, session = _here_and_other(tmp_path, monkeypatch)
    words = "words for the run's thread"
    held, answered = [], []
    page, errors, close = _page_on(session, [dict(_CONFIG)])
    try:
        page.on("response", lambda r: answered.append(r) if "/api/agent/run" in r.url else None)
        _say(page, "chignolin")
        page.wait_for_selector(_RUN)
        _shows(page, study_key(here))
        page.wait_for_timeout(800)
        page.route("**/api/status*", lambda route: held.append(route))
        page.click(_RUN)
        page.focus("#agent-request")
        page.keyboard.type(words)
        time.sleep(1.0)
        _switch_elsewhere(session, other)
        page.unroute("**/api/status*")
        for route in list(held):
            route.continue_()
        _shows(page, study_key(other))
        while not answered:
            page.wait_for_timeout(100)
        page.wait_for_timeout(2000)
        box = page.input_value("#agent-request")
        kept = _unsent(page)
        note = _note(page)
    finally:
        close()
        session.server.shutdown()
    run = study_key(workspace / "run1")
    assert made
    assert kept.get(run) == words, kept
    assert words not in kept.get(study_key(here), ""), kept
    assert words not in box
    assert "kept for it, not put here" in note, note
    assert errors == []


#: Every request of the dashboard's poll: held together, the page reads
#: no part of a poll without the others (CI-774 review L3).
_POLL = re.compile(r"/api/(app-state|status|metrics|events|results|structure-info|frames-info)(\?|$)")


def _hold(page, held):
    page.route(_POLL, lambda route: held.append(route))


def _release(page, held):
    page.unroute(_POLL)
    for route in list(held):
        route.continue_()
    del held[:]


def _texts(workspace, part):
    kept = _entries_everywhere(workspace)
    return [e.get("text") for path, entries in kept.items() if part in path for e in entries]


def _asking(page, asked):
    page.on("request", lambda r: asked.append(json.loads(r.post_data or "{}"))
            if "/api/agent/propose" in r.url else None)


def _events(page):
    page.evaluate("""() => { window.__else = [];
        window.addEventListener('dashboard:study-changed-elsewhere',
                                (e) => window.__else.push(e.detail.to)); }""")


#: Said beside the box while Send waits (CI-774 round 3: a message held to
#: be sent later by itself was lost, sent twice or sent elsewhere).
_SEND_WAITS = "Send waits until the run has started and this page shows which study is open."

_OFFERED = "() => document.getElementById('agent-propose').getAttribute('aria-disabled') === 'false'"


def _send_offered(page):
    return page.evaluate(_OFFERED)


def _send_when_offered(page, timeout=30000):
    """Once the page knows where a message would be answered, Send is
    offered again, and the person presses Enter on what is in the box (if
    the Enter before sent it, nothing is left to send)."""
    page.wait_for_function("() => (" + _OFFERED + ")()"
                           " || !document.getElementById('agent-request').value.trim()",
                           timeout=timeout)
    if page.input_value("#agent-request").strip():
        page.focus("#agent-request")
        page.keyboard.press("Enter")


def _requests(asked):
    return [a.get("request") for a in asked]


def _typed_and_entered_during_a_launch(page, here, words):
    """"chignolin" answered with a config; Run pressed; `words` typed and
    Enter pressed while the launch goes on: Send waits, and says so."""
    from fastmdxplora.gui.workspace import study_key

    _say(page, "chignolin")
    page.wait_for_selector(_RUN)
    _shows(page, study_key(here))
    page.wait_for_timeout(800)
    page.click(_RUN)
    page.focus("#agent-request")
    page.keyboard.type(words)
    page.keyboard.press("Enter")
    page.wait_for_timeout(300)
    return {"note": _note(page), "offered": _send_offered(page),
            "box": page.input_value("#agent-request")}


def test_a_run_s_study_passed_while_launching_is_not_this_window_s(tmp_path, monkeypatch):
    """L3 2, order A: the GUI goes to run1, then another window opens
    `other`, seen while the launch is unanswered. After the answer the page
    reads `other` again; the person types for it; another window opens run1.
    It is said as opened elsewhere, and the words stay kept for `other`."""
    from fastmdxplora.gui.workspace import study_key

    workspace, here, other, made, session = _here_and_other(tmp_path, monkeypatch, taking=4.0)
    run = workspace / "run1"
    held, answered = [], []
    page, errors, close = _page_on(session, [dict(_CONFIG)])
    try:
        page.on("response", lambda r: answered.append(r) if "/api/agent/run" in r.url else None)
        _say(page, "chignolin")
        page.wait_for_selector(_RUN)
        _shows(page, study_key(here))
        page.wait_for_timeout(800)
        _events(page)
        _hold(page, held)
        page.click(_RUN)
        time.sleep(1.0)
        _switch_elsewhere(session, other)
        _release(page, held)
        _shows(page, study_key(other))
        assert not answered
        while not answered:
            page.wait_for_timeout(100)
        page.wait_for_timeout(4000)
        page.fill("#agent-request", "typed for other")
        page.wait_for_timeout(600)
        _switch_elsewhere(session, run)
        _shows(page, study_key(run))
        page.wait_for_timeout(3000)
        said = page.evaluate("() => window.__else")
        box = page.input_value("#agent-request")
        kept = _unsent(page)
    finally:
        close()
        session.server.shutdown()
    assert made
    assert said == [study_key(other), study_key(run)], said
    assert "typed for other" not in box
    assert kept.get(study_key(other)) == "typed for other", kept
    assert errors == []


def test_an_agent_launch_answered_after_a_builder_launch_leaves_no_run_ahead(tmp_path,
                                                                           monkeypatch):
    """L3 2, order B: the Agent's Run (GUI to run1, answered after 4 s),
    then the Config Builder's Start (GUI to run2, answered at once). The
    page sees only run2. Another window then opens run1: said as opened
    elsewhere."""
    from fastmdxplora.gui.exploration import DashboardRuntime
    from fastmdxplora.gui.workspace import study_key

    workspace, here, other, made, session = _here_and_other(tmp_path, monkeypatch, taking=0.0)
    takes = [4.0, 0.3]

    def launch(self, state, *, dashboard_url=None, config=None):
        folder = workspace / f"run{len(made) + 1}"
        folder.mkdir()
        (folder / "exploration.yml").write_text("systems: []\n", encoding="utf-8")
        made.append(str(folder))
        self.active_root = folder
        time.sleep(takes[len(made) - 1])
        return {"ok": True, "output": str(folder)}

    monkeypatch.setattr(DashboardRuntime, "launch_from_config", launch)
    held, answered, built = [], [], []
    page, errors, close = _page_on(session, [dict(_CONFIG)])
    try:
        page.on("response", lambda r: answered.append(r) if "/api/agent/run" in r.url else None)
        page.on("response", lambda r: built.append(r) if r.url.endswith("/api/run") else None)
        _say(page, "chignolin")
        page.wait_for_selector(_RUN)
        _shows(page, study_key(here))
        page.wait_for_timeout(800)
        _events(page)
        _hold(page, held)
        page.click(_RUN)
        while not made:
            time.sleep(0.05)
        _builder(page, {"root": workspace})
        page.click("#run-start-button")
        while not (built and answered):
            page.wait_for_timeout(100)
        _release(page, held)
        run1, run2 = workspace / "run1", workspace / "run2"
        _shows(page, study_key(run2))
        page.wait_for_timeout(2000)
        _switch_elsewhere(session, run1)
        _shows(page, study_key(run1))
        page.wait_for_timeout(3000)
        said = page.evaluate("() => window.__else")
    finally:
        close()
        session.server.shutdown()
    assert len(made) == 2
    assert said == [study_key(run1)], said
    assert errors == []


def test_a_study_opened_and_gone_past_as_the_run_starts_is_said(tmp_path, monkeypatch):
    """L3 5: Run pressed; before the run's switch, another window opens
    `other`; the page sees `other`, then the run's study, both while the
    launch is unanswered. `other`'s opening is said; the run's is not."""
    from fastmdxplora.gui.workspace import study_key

    workspace, here, other, made, session = _here_and_other(tmp_path, monkeypatch, taking=6.0)
    runs, answered = [], []
    page, errors, close = _page_on(session, [dict(_CONFIG)])
    try:
        page.on("response", lambda r: answered.append(r) if "/api/agent/run" in r.url else None)
        _say(page, "chignolin")
        page.wait_for_selector(_RUN)
        _shows(page, study_key(here))
        page.wait_for_timeout(800)
        _events(page)
        page.evaluate(_HEARD)
        page.route("**/api/agent/run", lambda route: runs.append(route))
        page.click(_RUN)
        while not runs:
            page.wait_for_timeout(50)
        _switch_elsewhere(session, other)
        _shows(page, study_key(other))
        page.wait_for_timeout(500)
        runs[0].continue_()
        page.unroute("**/api/agent/run")
        _shows(page, study_key(workspace / "run1"))
        assert not answered
        while not answered:
            page.wait_for_timeout(100)
        page.wait_for_timeout(3000)
        said = page.evaluate("() => window.__else")
        heard = page.evaluate("() => window.__said")
    finally:
        close()
        session.server.shutdown()
    assert made
    assert said == [study_key(other)], said
    assert [h for h in heard if "was opened elsewhere as the run started" in h], heard
    assert errors == []


def test_words_kept_for_a_study_are_put_back_when_it_is_opened_again(tmp_path, monkeypatch):
    """L5 1 and L6 H1, with no launch: words typed on `here`, another window
    opens `other` (the words kept for `here`); words typed on `other`,
    another window opens `here` again. `here`'s words are back in its box,
    and `other`'s are kept for it: they were erased as the box emptied. The
    two are of one system, and the note names the study with its folder
    (L6 L1, L7 6)."""
    from fastmdxplora.gui.workspace import study_key

    workspace, here, other, made, session = _here_and_other(tmp_path, monkeypatch)
    for study in (here, other):
        (study / "manifest.json").write_text(json.dumps({"system": "tri-ala.pdb"}),
                                             encoding="utf-8")
    page, errors, close = _page_on(session, [])
    try:
        _shows(page, study_key(here))
        page.wait_for_timeout(800)
        page.fill("#agent-request", "for here")
        page.wait_for_timeout(500)
        _switch_elsewhere(session, other)
        _shows(page, study_key(other))
        page.wait_for_timeout(1000)
        note = _note(page)
        page.fill("#agent-request", "for other")
        page.wait_for_timeout(500)
        _switch_elsewhere(session, here)
        _shows(page, study_key(here))
        page.wait_for_timeout(1500)
        box = page.input_value("#agent-request")
        kept = _unsent(page)
    finally:
        close()
        session.server.shutdown()
    assert box == "for here", (box, kept)
    assert kept.get(study_key(other)) == "for other", kept
    assert "is kept for it" in note and "(here)" in note, note
    assert errors == []


def test_the_read_asked_at_a_launch_s_answer_counts_as_after_it(tmp_path, monkeypatch):
    """L5 3 and L6 M4: the status read asked as a launch answers is asked in
    the same millisecond; counted in polls, it is one asked after the
    answer, and the study is known once it is applied."""
    from fastmdxplora.gui.workspace import study_key

    workspace, here, other, made, session = _here_and_other(tmp_path, monkeypatch)
    page, errors, close = _page_on(session, [])
    try:
        _shows(page, study_key(here))
        page.wait_for_timeout(1500)
        known = []
        for _ in range(3):
            page.evaluate("""() => { const d = window.FastMDXDashboard;
                d.studyChangeAsked(); d.studyChangeAnswered(""); }""")
            page.wait_for_timeout(1500)
            known.append(page.evaluate("() => window.FastMDXDashboard.studyKnown()"))
    finally:
        close()
        session.server.shutdown()
    assert known == [True, True, True], known
    assert errors == []


def _stale_status(page, stale):
    """The next status read is answered as the server reads it now, and
    handed to the page only when the test says (`stale[0][0].fulfill`)."""
    def hold(route):
        if stale:
            route.continue_()
        else:
            stale.append((route, route.fetch()))

    page.route("**/api/status*", hold)
    page.evaluate("() => document.getElementById('refresh-now').click()")
    while not stale:
        page.wait_for_timeout(50)


def test_a_read_of_another_study_asked_before_the_answer_leaves_the_run_this_window_s(
        tmp_path, monkeypatch):
    """L6 M3: Run pressed; another window opens `other` before the run's
    switch; the page's read of `other`, asked then, is applied after the
    answer. `other` is said as opened elsewhere; the run's own study, seen
    next, is not."""
    from fastmdxplora.gui.workspace import study_key

    workspace, here, other, made, session = _here_and_other(tmp_path, monkeypatch, taking=0.0)
    runs, stale, answered = [], [], []
    page, errors, close = _page_on(session, [dict(_CONFIG)])
    try:
        page.on("response", lambda r: answered.append(r) if "/api/agent/run" in r.url else None)
        _say(page, "chignolin")
        page.wait_for_selector(_RUN)
        _shows(page, study_key(here))
        page.wait_for_timeout(800)
        _events(page)
        page.route("**/api/agent/run", lambda route: runs.append(route))
        page.click(_RUN)
        while not runs:
            page.wait_for_timeout(50)
        _switch_elsewhere(session, other)
        _stale_status(page, stale)
        runs[0].continue_()
        page.unroute("**/api/agent/run")
        while not answered:
            page.wait_for_timeout(50)
        page.wait_for_timeout(200)
        stale[0][0].fulfill(response=stale[0][1])
        _shows(page, study_key(workspace / "run1"))
        page.wait_for_timeout(3000)
        said = page.evaluate("() => window.__else")
    finally:
        close()
        session.server.shutdown()
    assert made
    assert said == [study_key(other)], said
    assert errors == []


def test_an_unreadable_answer_takes_the_run_s_study_seen_last_as_its_own(tmp_path, monkeypatch):
    """L5 4: Run pressed; another window opens `other`, which the page sees;
    the run starts (the GUI goes to run1, which the page sees) and its answer
    cannot be read. `other` is said as opened elsewhere; run1 is not."""
    from fastmdxplora.gui.workspace import study_key

    workspace, here, other, made, session = _here_and_other(tmp_path, monkeypatch, taking=0.0)
    runs = []
    page, errors, close = _page_on(session, [dict(_CONFIG)])
    try:
        _say(page, "chignolin")
        page.wait_for_selector(_RUN)
        _shows(page, study_key(here))
        page.wait_for_timeout(800)
        _events(page)
        page.route("**/api/agent/run", lambda route: runs.append(route))
        page.click(_RUN)
        while not runs:
            page.wait_for_timeout(50)
        _switch_elsewhere(session, other)
        _shows(page, study_key(other))
        runs[0].fetch()
        _shows(page, study_key(workspace / "run1"))
        runs[0].fulfill(status=502, content_type="text/html", body="<h1>Bad gateway</h1>")
        page.wait_for_selector("text=I could not read the software's answer", timeout=30000)
        page.wait_for_timeout(2000)
        said = page.evaluate("() => window.__else")
    finally:
        close()
        session.server.shutdown()
    assert made
    assert said == [study_key(other)], said
    assert errors == []


def test_enter_during_a_launch_sends_nothing_until_the_person_sends_again(tmp_path,
                                                                           monkeypatch):
    """L10 F1 and round 3: Enter during the launch sends nothing; it says
    why beside the box, and the words stay there. They go with the thread
    to the run's study; once the page shows it, Send is offered again, and
    nothing is sent until the person sends it, into the run's thread with
    the run's view. A message held to be sent later by itself was lost,
    sent twice, or sent into another thread or study."""
    from fastmdxplora.gui.workspace import study_key

    workspace, here, other, made, session = _here_and_other(tmp_path, monkeypatch, taking=3.0)
    run = workspace / "run1"
    words = "typed as it starts"
    asked, answered = [], []
    page, errors, close = _page_on(session, [dict(_CONFIG), _answer("Answered.")])
    try:
        _asking(page, asked)
        page.on("response", lambda r: answered.append(r) if "/api/agent/run" in r.url else None)
        during = _typed_and_entered_during_a_launch(page, here, words)
        while not answered:
            page.wait_for_timeout(50)
        _shows(page, study_key(run))
        page.wait_for_function(_OFFERED, timeout=15000)
        page.wait_for_timeout(3000)
        before, note = _requests(asked), _note(page)
        _send_when_offered(page)
        page.wait_for_selector("#agent-thread .agent-answer:has-text('Answered.')", timeout=15000)
    finally:
        close()
        session.server.shutdown()
    assert during == {"note": _SEND_WAITS, "offered": False, "box": words}, during
    assert before == ["chignolin"] and note == "", (before, note)
    assert _requests(asked) == ["chignolin", words]
    assert (asked[1].get("current_view") or {}).get("study", "").endswith("run1"), asked[1]
    assert words in _texts(workspace, "run1/agent")
    assert errors == []


def test_a_study_seen_first_takes_nothing_sent_and_blocks_nothing(tmp_path, monkeypatch):
    """L3 1, L4 H1, L5 2 and L6 M2: Enter during the launch; the launch
    answers before the page reads a status, and another window opens
    `other`, which the page sees first. Nothing is sent, into the run's
    thread or `other`'s; the words are kept for the run's study, where the
    thread they were typed in went. On `other`, a question is sent at once."""
    from fastmdxplora.gui.workspace import study_key

    workspace, here, other, made, session = _here_and_other(tmp_path, monkeypatch, taking=3.0)
    words = "typed as it starts"
    held, answered, asked = [], [], []
    page, errors, close = _page_on(session, [dict(_CONFIG), _answer("About other.")])
    try:
        page.on("response", lambda r: answered.append(r) if "/api/agent/run" in r.url else None)
        _asking(page, asked)
        _hold(page, held)
        during = _typed_and_entered_during_a_launch(page, here, words)
        while not answered:
            page.wait_for_timeout(50)
        page.wait_for_timeout(400)
        _switch_elsewhere(session, other)
        _release(page, held)
        _shows(page, study_key(other))
        page.wait_for_timeout(3000)
        box, kept, before = page.input_value("#agent-request"), _unsent(page), _requests(asked)
        page.fill("#agent-request", "a question about other")
        page.keyboard.press("Enter")
        page.wait_for_selector("#agent-thread .agent-answer:has-text('About other.')", timeout=15000)
    finally:
        close()
        session.server.shutdown()
    assert during["note"] == _SEND_WAITS, during
    assert before == ["chignolin"], before
    assert words not in box and kept.get(study_key(workspace / "run1")) == words, (box, kept)
    assert _requests(asked) == ["chignolin", "a question about other"]
    assert (asked[1].get("current_view") or {}).get("study", "").endswith("other")
    assert words not in _texts(workspace, "run1/agent") + _texts(workspace, "other/agent")
    assert errors == []


def test_words_that_follow_the_thread_are_sent_only_by_the_person(tmp_path, monkeypatch):
    """L4 M1, L3 6 and L5 5, in L4's order: Enter during the launch; another
    window opens `other` before the run's switch and the page sees it; the
    run starts, and the page sees its study. The words went with the thread
    there; Send is offered soon after, and they are sent when the person
    sends them, with the run's view."""
    from fastmdxplora.gui.workspace import study_key

    workspace, here, other, made, session = _here_and_other(tmp_path, monkeypatch, taking=0.3)
    run = workspace / "run1"
    words = "words for the run's thread"
    runs, asked = [], []
    page, errors, close = _page_on(session, [dict(_CONFIG), _answer("Answered.")])
    try:
        _asking(page, asked)
        page.route("**/api/agent/run", lambda route: runs.append(route))
        during = _typed_and_entered_during_a_launch(page, here, words)
        while not runs:
            page.wait_for_timeout(50)
        _switch_elsewhere(session, other)
        _shows(page, study_key(other))
        page.wait_for_timeout(500)
        runs[0].continue_()
        page.unroute("**/api/agent/run")
        _shows(page, study_key(run))
        page.wait_for_function(_OFFERED, timeout=6000)
        page.wait_for_timeout(2000)
        before, box = _requests(asked), page.input_value("#agent-request")
        _send_when_offered(page)
        page.wait_for_selector("#agent-thread .agent-answer:has-text('Answered.')", timeout=15000)
        note = _note(page)
    finally:
        close()
        session.server.shutdown()
    assert during["note"] == _SEND_WAITS, during
    assert before == ["chignolin"] and box == words, (before, box)
    assert (asked[-1].get("current_view") or {}).get("study", "").endswith("run1")
    assert words in _texts(workspace, "run1/agent") and note == "", note
    assert errors == []


def test_an_unreadable_answer_sends_nothing_by_itself(tmp_path, monkeypatch):
    """L3 3: Enter during the launch; the server starts the run (the GUI
    goes to run1, the thread with it) but its answer cannot be read.
    Nothing is sent by itself; once the page has read which study is open,
    the person sends it, with that study open."""
    import urllib.request

    from fastmdxplora.gui.workspace import study_key

    workspace, here, other, made, session = _here_and_other(tmp_path, monkeypatch, taking=1.0)
    words = "typed as it starts"
    held, asked = [], []
    base = session.url.split("#")[0].rstrip("/")

    def gui_study():
        with urllib.request.urlopen(base + "/api/status", timeout=10) as answer:
            return json.loads(answer.read()).get("study")

    def lost(route):
        route.fetch()
        route.fulfill(status=502, content_type="text/html", body="<h1>Bad gateway</h1>")

    page, errors, close = _page_on(session, [dict(_CONFIG), _answer("Answered.")])
    try:
        page.on("request", lambda r: asked.append((json.loads(r.post_data or "{}"), gui_study()))
                if "/api/agent/propose" in r.url else None)
        page.route("**/api/agent/run", lost)
        _hold(page, held)
        during = _typed_and_entered_during_a_launch(page, here, words)
        page.wait_for_selector("text=I could not read the software's answer", timeout=30000)
        page.wait_for_timeout(1500)
        _release(page, held)
        _shows(page, study_key(workspace / "run1"))
        page.wait_for_function(_OFFERED, timeout=15000)
        page.wait_for_timeout(2000)
        before = [a.get("request") for a, _ in asked]
        _send_when_offered(page)
        page.wait_for_selector("#agent-thread .agent-answer:has-text('Answered.')", timeout=30000)
    finally:
        close()
        session.server.shutdown()
    assert during["note"] == _SEND_WAITS, during
    assert before == ["chignolin"], before
    sent = [(a, gui) for a, gui in asked if a.get("request") == words]
    assert len(sent) == 1, asked
    view = (sent[0][0].get("current_view") or {}).get("study") or ""
    assert study_key(view) == sent[0][1] == study_key(workspace / "run1"), (view, sent[0][1])
    assert words in _texts(workspace, "run1/agent")
    assert errors == []


def test_a_stale_read_of_the_start_study_leaves_the_words_with_the_run(tmp_path, monkeypatch):
    """L6 M1: Enter during the launch; its answer cannot be read; a status
    read of the start study, asked before the run's switch, is applied
    after the answer. The words still go with the thread to the run's
    study, and are sent there when the person sends them."""
    from fastmdxplora.gui.workspace import study_key

    workspace, here, other, made, session = _here_and_other(tmp_path, monkeypatch, taking=1.0)
    words = "typed as it starts"
    stale, asked = [], []

    def lost(route):
        route.fetch()
        route.fulfill(status=502, content_type="text/html", body="<h1>Bad gateway</h1>")

    page, errors, close = _page_on(session, [dict(_CONFIG), _answer("Answered.")])
    try:
        _asking(page, asked)
        _say(page, "chignolin")
        page.wait_for_selector(_RUN)
        _shows(page, study_key(here))
        page.wait_for_timeout(800)
        page.route("**/api/agent/run", lost)
        _stale_status(page, stale)
        page.click(_RUN)
        page.focus("#agent-request")
        page.keyboard.type(words)
        page.keyboard.press("Enter")
        page.wait_for_selector("text=I could not read the software's answer", timeout=30000)
        stale[0][0].fulfill(response=stale[0][1])
        _shows(page, study_key(workspace / "run1"))
        _send_when_offered(page)
        page.wait_for_selector("#agent-thread .agent-answer:has-text('Answered.')", timeout=30000)
        kept = _unsent(page)
    finally:
        close()
        session.server.shutdown()
    sent = [a for a in asked if a.get("request") == words]
    assert len(sent) == 1 and (sent[0].get("current_view") or {}).get("study", "").endswith("run1")
    assert words in _texts(workspace, "run1/agent")
    assert words not in json.dumps(kept), kept
    assert errors == []


def test_a_chat_of_no_study_waits_and_is_filed_as_one(tmp_path, monkeypatch):
    """L3 4 and L9 2: after a launch answers, before the page has read which
    study is open, the person starts a chat of no study and presses Enter:
    nothing is sent, as for any message then. Sent by the person once Send
    is offered, it is the chat's, never filed under a study."""
    from fastmdxplora.gui.workspace import study_key

    workspace, here, other, made, session = _here_and_other(tmp_path, monkeypatch, taking=3.0)
    words = "a general question"
    held, answered, asked = [], [], []
    page, errors, close = _page_on(session, [dict(_CONFIG), _answer("General answer.")])
    try:
        page.on("response", lambda r: answered.append(r) if "/api/agent/run" in r.url else None)
        _asking(page, asked)
        _say(page, "chignolin")
        page.wait_for_selector(_RUN)
        _shows(page, study_key(here))
        page.wait_for_timeout(800)
        _hold(page, held)
        page.click(_RUN)
        time.sleep(1.0)
        _switch_elsewhere(session, here)
        while not answered:
            page.wait_for_timeout(50)
        page.wait_for_timeout(300)
        page.evaluate("window.FastMDXAgentPanel.fresh(null)")
        page.wait_for_timeout(800)
        page.fill("#agent-request", words)
        page.keyboard.press("Enter")
        page.wait_for_timeout(500)
        note, before = _note(page), _requests(asked)
        _release(page, held)
        _send_when_offered(page)
        page.wait_for_selector("#agent-thread .agent-answer:has-text('General answer.')",
                               timeout=15000)
        page.wait_for_timeout(1000)
    finally:
        close()
        session.server.shutdown()
    assert note == _SEND_WAITS and before == ["chignolin"], (note, before)
    assert _requests(asked) == ["chignolin", words]
    assert words not in _texts(workspace, "here/agent") + _texts(workspace, "run1/agent")
    assert words in _texts(workspace, "")
    assert errors == []


def test_a_message_waiting_is_never_sent_as_the_page_goes(tmp_path, monkeypatch):
    """L4 H2 and L8 H1: Enter during the launch, then the page goes (its
    pagehide; the page is kept to watch what follows). Nothing is sent as
    it goes, and the words are kept as any unsent words."""
    from fastmdxplora.gui.workspace import study_key

    workspace, here, other, made, session = _here_and_other(tmp_path, monkeypatch, taking=2.0)
    words = "typed, then the page went"
    asked, answered = [], []
    page, errors, close = _page_on(session, [dict(_CONFIG), _answer("Answered.")])
    try:
        page.on("response", lambda r: answered.append(r) if "/api/agent/run" in r.url else None)
        _asking(page, asked)
        _typed_and_entered_during_a_launch(page, here, words)
        page.evaluate("() => window.dispatchEvent(new PageTransitionEvent('pagehide', {persisted: false}))")
        while not answered:
            page.wait_for_timeout(50)
        _shows(page, study_key(workspace / "run1"))
        page.evaluate("() => window.FastMDXDashboard.readStatusNow()")
        page.wait_for_timeout(3000)
        kept = _unsent(page)
    finally:
        close()
        session.server.shutdown()
    assert made
    assert _requests(asked) == ["chignolin"], asked
    assert words in json.dumps(kept), kept
    assert errors == []


def test_reloading_during_a_launch_asks_nothing_and_sends_nothing(tmp_path, monkeypatch):
    """L7 2, L7 4, L8 H1 and L10 F1: Enter during the launch, then a reload.
    The words are kept in this tab as any unsent words, so the browser
    does not ask; nothing is ever sent by itself after the reload."""
    from fastmdxplora.gui.workspace import study_key

    workspace, here, other, made, session = _here_and_other(tmp_path, monkeypatch, taking=3.0)
    words = "typed, then reloaded"
    asked, dialogs = [], []
    page, errors, close = _page_on(session, [dict(_CONFIG), _answer("Answered.")])
    try:
        _asking(page, asked)
        page.on("dialog", lambda d: (dialogs.append(d.type), d.accept()))
        _typed_and_entered_during_a_launch(page, here, words)
        page.reload(wait_until="domcontentloaded")
        _shows(page, study_key(workspace / "run1"))
        page.wait_for_timeout(5000)
        kept = _unsent(page)
        box = page.input_value("#agent-request")
    finally:
        close()
        session.server.shutdown()
    assert dialogs == [], dialogs
    assert _requests(asked) == ["chignolin"], asked
    assert words in json.dumps(kept) or words in box, (kept, box)
    assert errors == []


def test_the_page_s_own_switches_go_on_with_words_kept(tmp_path, monkeypatch):
    """L9 1: words kept for the run's study (Enter during the launch,
    another window opening `other` first). The page opens `here`, then the
    run's study, by their tabs: no "Leave site?", each switch made; the
    words are in the run's box, and nothing is sent."""
    from fastmdxplora.gui.workspace import study_key

    workspace, here, other, made, session = _here_and_other(tmp_path, monkeypatch, taking=3.0)
    run = workspace / "run1"
    words = "typed as it starts"
    held, answered, asked, dialogs = [], [], [], []
    page, errors, close = _page_on(session, [dict(_CONFIG)])
    try:
        page.on("response", lambda r: answered.append(r) if "/api/agent/run" in r.url else None)
        page.on("dialog", lambda d: (dialogs.append(d.type), d.dismiss()))
        _asking(page, asked)
        _hold(page, held)
        _typed_and_entered_during_a_launch(page, here, words)
        while not answered:
            page.wait_for_timeout(50)
        page.wait_for_timeout(400)
        _switch_elsewhere(session, other)
        _release(page, held)
        _shows(page, study_key(other))
        page.wait_for_timeout(2000)
        for study in (here, run):
            page.evaluate("(key) => { window.FastMDXDashboard.openStudyTab(key, 'agent'); }",
                          study_key(study))
            page.wait_for_timeout(1500)
            page.wait_for_function("() => window.FastMDXDashboard && window.FastMDXDashboard.state")
            _shows(page, study_key(study))
        page.wait_for_selector("#agent-request", state="visible")
        page.wait_for_timeout(1500)
        box = page.input_value("#agent-request")
    finally:
        close()
        session.server.shutdown()
    assert dialogs == [], dialogs
    assert box == words, box
    assert _requests(asked) == ["chignolin"], asked
    assert errors == []


def test_a_retry_during_a_launch_sends_nothing(tmp_path, monkeypatch):
    """L9 3 and L10 F4: "hello" answered; Run pressed; Retry on "hello"
    while the launch goes on: nothing is sent, now or later by itself, and
    Send waits, saying so."""
    from tests.test_a_conversation_follows_a_run_that_wrote_nothing import _RETRY
    from fastmdxplora.gui.workspace import study_key

    workspace, here, other, made, session = _here_and_other(tmp_path, monkeypatch, taking=3.0)
    asked, answered = [], []
    page, errors, close = _page_on(session, [_answer("Hello."), dict(_CONFIG), _answer("Again.")])
    try:
        _asking(page, asked)
        page.on("response", lambda r: answered.append(r) if "/api/agent/run" in r.url else None)
        _shows(page, study_key(here))
        _say(page, "hello")
        page.wait_for_selector("#agent-thread .agent-answer:has-text('Hello.')")
        _say(page, "chignolin")
        page.wait_for_selector(_RUN)
        page.wait_for_timeout(800)
        page.click(_RUN)
        hello = page.locator("#agent-thread .agent-msg-user").first
        hello.hover()
        hello.locator(_RETRY).first.click(force=True)
        page.wait_for_timeout(300)
        note = _note(page)
        while not answered:
            page.wait_for_timeout(50)
        _shows(page, study_key(workspace / "run1"))
        page.wait_for_timeout(4000)
    finally:
        close()
        session.server.shutdown()
    assert note == _SEND_WAITS, note
    assert _requests(asked) == ["hello", "chignolin"], asked
    assert errors == []


def test_a_config_builder_launch_makes_send_wait_too(tmp_path, monkeypatch):
    """L10 F3: while a launch of the Config Builder is unanswered (asked
    and answered as the Builder asks them), Enter sends nothing and says
    so; once it answers, Send is offered, and the message goes into the
    study shown's thread."""
    from fastmdxplora.gui.workspace import study_key

    workspace, here, other, made, session = _here_and_other(tmp_path, monkeypatch)
    words = "asked as the Builder starts"
    asked = []
    page, errors, close = _page_on(session, [_answer("Answered.")])
    try:
        _asking(page, asked)
        _shows(page, study_key(here))
        page.wait_for_timeout(1500)
        page.evaluate("() => window.FastMDXDashboard.studyChangeAsked()")
        page.fill("#agent-request", words)
        page.keyboard.press("Enter")
        page.wait_for_timeout(1000)
        note, offered, before = _note(page), _send_offered(page), _requests(asked)
        page.evaluate("() => window.FastMDXDashboard.studyChangeAnswered('')")
        _send_when_offered(page)
        page.wait_for_selector("#agent-thread .agent-answer:has-text('Answered.')", timeout=15000)
    finally:
        close()
        session.server.shutdown()
    assert (note, offered, before) == (_SEND_WAITS, False, []), (note, offered, before)
    assert _requests(asked) == [words]
    assert (asked[0].get("current_view") or {}).get("study", "").endswith("here")
    assert words in _texts(workspace, "here/agent")
    assert errors == []


def test_a_send_goes_with_the_view_read_after_its_study_s_status(tmp_path, monkeypatch):
    """L8 M2: one poll reads the app state with `here` open, another window
    opens `other`, and the poll's status is read after it: the page shows
    `other` with `here` as its view, and the notice named `here`. Enter
    sends once the view is read again after the status, with `other` as
    its view, and the notice names `other`."""
    from fastmdxplora.gui.workspace import study_key

    workspace, here, other, made, session = _here_and_other(tmp_path, monkeypatch)
    words = "about the study shown"
    asked, apps, later, statuses = [], [], [], []

    def app(route):
        if apps:
            later.append(route)
        else:
            apps.append((route, route.fetch()))

    def status(route):
        if statuses:
            route.continue_()
        else:
            statuses.append(route)

    page, errors, close = _page_on(session, [_answer("Answered.")])
    try:
        _asking(page, asked)
        _shows(page, study_key(here))
        page.wait_for_timeout(1500)
        page.evaluate(_HEARD)
        page.route("**/api/app-state*", app)
        page.route("**/api/status*", status)
        page.evaluate("() => document.getElementById('refresh-now').click()")
        while not (apps and statuses):
            page.wait_for_timeout(50)
        _switch_elsewhere(session, other)
        apps[0][0].fulfill(response=apps[0][1])
        statuses[0].continue_()
        page.unroute("**/api/status*")
        _shows(page, study_key(other))
        page.wait_for_timeout(500)
        page.fill("#agent-request", words)
        page.keyboard.press("Enter")
        page.wait_for_timeout(500)
        note, before = _note(page), _requests(asked)
        page.unroute("**/api/app-state*")
        for route in later:
            route.continue_()
        _send_when_offered(page)
        page.wait_for_selector("#agent-thread .agent-answer:has-text('Answered.')", timeout=15000)
        page.wait_for_function("() => window.__said.some((s) => /opened elsewhere/.test(s))",
                               timeout=15000)
        said = [s for s in page.evaluate("() => window.__said") if "opened elsewhere" in s]
    finally:
        close()
        session.server.shutdown()
    assert before == [], (note, before)
    assert len(asked) == 1, asked
    assert (asked[0].get("current_view") or {}).get("study", "").endswith("other"), asked
    assert said and said[-1].startswith("other was opened elsewhere"), said
    assert errors == []


def test_send_and_the_notice_do_not_wait_for_the_study_s_results(tmp_path, monkeypatch):
    """L10 F6: another window opens `other`, whose results take seconds to
    read. Once the page shows `other`, a message sent there goes as soon as
    the view is read after its status, and the notice is said: the study's
    results are not waited for again (they were read only for its name)."""
    from fastmdxplora.gui.workspace import study_key

    workspace, here, other, made, session = _here_and_other(tmp_path, monkeypatch)
    pending, asked = [], []
    page, errors, close = _page_on(session, [_answer("Answered.")])
    try:
        page.on("request", lambda r: asked.append(time.time())
                if "/api/agent/propose" in r.url else None)
        _shows(page, study_key(here))
        page.wait_for_timeout(1500)
        page.evaluate(_HEARD)
        page.route("**/api/results*", lambda route: pending.append((time.time(), route)))

        def pump(until, most=40.0):
            end = time.time() + most
            while time.time() < end:
                now = time.time()
                for when, route in list(pending):
                    if now - when > 5.0:
                        pending.remove((when, route))
                        route.continue_()
                if until():
                    return True
                page.wait_for_timeout(50)
            return False

        _switch_elsewhere(session, other)
        assert pump(lambda: page.evaluate("() => window.FastMDXDashboard.state.studyShown")
                    == study_key(other))
        shown = time.time()
        page.fill("#agent-request", "about other")
        page.keyboard.press("Enter")
        assert pump(lambda: bool(asked), 20.0)
        assert pump(lambda: page.evaluate(
            "() => window.__said.some((s) => /opened elsewhere/.test(s))"), 20.0)
        said = time.time()
        page.unroute("**/api/results*")
        for _, route in pending:
            route.continue_()
    finally:
        close()
        session.server.shutdown()
    assert asked[0] - shown < 3.0 and said - shown < 3.0, (asked[0] - shown, said - shown)
    assert errors == []


def test_a_config_builder_launch_takes_only_its_run_s_study_as_its_own(browser, site,  # noqa: F811
                                                                         monkeypatch):
    """L1 4 and L2 M2: the Config Builder's Start pressed, another window
    opens a study while the launch waits, and this page sees it; then the
    run starts and this page sees its study. The other window's study is
    said as opened elsewhere, the run's own is not."""
    from fastmdxplora.gui.workspace import study_key

    made = _launches(monkeypatch, site["root"])
    page = _open(browser, site)
    try:
        page.evaluate("""() => { window.__else = [];
            window.addEventListener('dashboard:study-changed-elsewhere',
                                    (e) => window.__else.push(e.detail.to)); }""")
        _builder(page, site)
        page.wait_for_timeout(800)
        held = []
        page.route("**/api/run", lambda route: held.append(route))
        with page.expect_request("**/api/run"):
            page.click("#run-start-button")
        _switch_elsewhere(site["session"], site["studies"][1])
        _shows(page, site["keys"][1])
        page.wait_for_timeout(500)
        held[0].continue_()
        page.unroute("**/api/run")
        _shows(page, study_key(site["root"] / "run1"))
        page.wait_for_timeout(3000)
        said = page.evaluate("() => window.__else")
        errors = list(page.errors)
    finally:
        page.context.close()
    assert made
    assert said == [site["keys"][1]], said
    assert errors == []


def test_run_as_it_is_takes_its_run_s_study_as_its_own(browser, site, monkeypatch):  # noqa: F811
    """L4 M2: the Config Builder's "Run on this machine as it is" opens its
    run's study as this window's own change: not said as opened elsewhere."""
    from fastmdxplora.gui.exploration import DashboardRuntime
    from fastmdxplora.gui.workspace import study_key

    root = site["root"]
    config = root / "as_it_is.yml"
    config.write_text("systems:\n- system: 1UAO\n", encoding="utf-8")
    made = []

    def launch(self, config_path, output=None, dashboard_url=None):
        folder = root / "run1"
        folder.mkdir()
        (folder / "exploration.yml").write_text("systems: []\n", encoding="utf-8")
        made.append(str(folder))
        self.active_root = folder
        return {"ok": True, "output": str(folder), "config_path": str(config_path)}

    monkeypatch.setattr(DashboardRuntime, "launch_existing_config", launch)
    page = _open(browser, site)
    try:
        page.route("**/api/check-config", lambda route: route.fulfill(
            status=200, content_type="application/json",
            body=json.dumps({"ok": True, "phases": ["simulation"], "systems": 1,
                             "settings_named": 1})))
        _events(page)
        page.evaluate("window.FastMDXDashboard.navigate('run')")
        page.evaluate("() => window.FastMDXRun.setStart('config')")
        page.fill("#run-config-path", str(config))
        page.click("#run-check-config")
        page.wait_for_function("() => !document.getElementById('run-as-is').disabled")
        page.click("#run-as-is")
        _shows(page, study_key(root / "run1"))
        page.wait_for_timeout(3000)
        said = page.evaluate("() => window.__else")
        errors = list(page.errors)
    finally:
        page.context.close()
    assert made
    assert said == [], said
    assert errors == []


def test_an_older_poll_answered_late_is_not_applied_over_a_newer_one(browser, site):  # noqa: F811
    """L4 L2: a poll's status read before another window opened a study,
    answered after a later poll that read that study: not applied, so the
    page neither goes back nor says the study it showed opened elsewhere."""
    page = _open(browser, site)
    late = []
    try:
        page.wait_for_timeout(1000)
        _events(page)

        def slow(route):
            if late:
                route.continue_()
            else:
                late.append((route, route.fetch()))

        page.route("**/api/status*", slow)
        page.evaluate("() => document.getElementById('refresh-now').click()")
        while not late:
            page.wait_for_timeout(50)
        _switch_elsewhere(site["session"], site["studies"][1])
        page.evaluate("() => document.getElementById('refresh-now').click()")
        _shows(page, site["keys"][1])
        route, answer = late[0]
        route.fulfill(response=answer)
        page.wait_for_timeout(2000)
        shown = page.evaluate("() => window.FastMDXDashboard.state.studyShown")
        said = page.evaluate("() => window.__else")
        errors = list(page.errors)
    finally:
        page.context.close()
    assert shown == site["keys"][1]
    assert said == [site["keys"][1]], said
    assert errors == []


def test_run_as_it_is_with_an_unreadable_answer_is_offered_again(browser, site):  # noqa: F811
    """L6 L2: "Run on this machine as it is" whose answer cannot be read:
    said as Start says it, its button offered again, no page error."""
    root = site["root"]
    config = root / "as_it_is.yml"
    config.write_text("systems:\n- system: 1UAO\n", encoding="utf-8")
    page = _open(browser, site)
    try:
        page.route("**/api/check-config", lambda route: route.fulfill(
            status=200, content_type="application/json",
            body=json.dumps({"ok": True, "phases": ["simulation"], "systems": 1,
                             "settings_named": 1})))
        page.route("**/api/run-config", lambda route: route.fulfill(
            status=502, content_type="text/html", body="<h1>Bad gateway</h1>"))
        page.evaluate("window.FastMDXDashboard.navigate('run')")
        page.evaluate("() => window.FastMDXRun.setStart('config')")
        page.fill("#run-config-path", str(config))
        page.click("#run-check-config")
        page.wait_for_function("() => !document.getElementById('run-as-is').disabled")
        page.click("#run-as-is")
        page.wait_for_function(
            "() => document.getElementById('run-note').textContent === 'Could not reach the server.'")
        page.wait_for_function("() => !document.getElementById('run-as-is').disabled")
        page.wait_for_timeout(500)
        errors = list(page.errors)
    finally:
        page.context.close()
    assert errors == []
