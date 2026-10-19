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
    """K5 1, with the send held: run A answered and not yet seen, a message
    asking for run B waits until the page has read which study is open (it
    sees A); B is started, and the words typed go with the thread to B, and
    are kept for B."""
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
        # The page does not see run A (run2) yet.
        page.route("**/api/status*", lambda route: statuses.append(route))
        page.click(f"{_RUN}:visible >> nth=-1")
        page.wait_for_function(
            "() => document.querySelectorAll('#agent-thread .agent-running').length >= 2",
            timeout=30000)
        _say(page, "and again")
        page.wait_for_timeout(1500)
        # Held: the thread went to run A, which the page has not seen.
        held = _note(page)
        cards = page.locator("#agent-thread .agent-study:not([hidden])").count()
        assert page.evaluate("() => window.FastMDXDashboard.state.studyShown") == run1
        page.unroute("**/api/status*")
        for route in list(statuses):
            route.continue_()
        _shows(page, run2)
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
    assert cards == 2 and "Sent once the run has started" in held, (cards, held)
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
    it answers, before the page's next status. It is answered with the
    study shown open and kept in that study's thread, never in the run's."""
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
        if sent != "during the launch":
            page.wait_for_timeout(1000)
            page.unroute("**/api/status*")
            for route in held[len(released):]:
                route.continue_()
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


def test_a_launched_draft_is_not_written_back_on_reload(browser, site):  # noqa: F811
    """L1 3 and L2 L2: a study started from the Config Builder leaves no
    draft of itself: leaving the page wrote the launched form back."""
    page = _open(browser, site)
    try:
        _builder(page, site)
        page.wait_for_timeout(800)
        out = str(site["root"] / "launched")
        page.route("**/api/run", lambda route: route.fulfill(
            status=200, content_type="application/json",
            body=json.dumps({"ok": True, "output": out, "config_path": out + "/config.yml"})))
        page.click("#run-start-button")
        page.wait_for_function("() => /Running/.test(document.getElementById('run-note').textContent)")
        page.reload()
        page.wait_for_function("() => window.FastMDXRun && window.FastMDXRun.state.schema")
        page.wait_for_timeout(1000)
        drafts = page.evaluate("""() => Object.keys(localStorage)
            .filter((key) => key.startsWith('fmx.builderDraft'))""")
    finally:
        page.context.close()
    assert drafts == []
