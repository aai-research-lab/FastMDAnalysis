"""A file opened in the reader as the page loads is said as of its own study.

The reader keeps the study a file is of, so that a file read in one study is
not said to the Agent as open in another. It took that study from the page's
app state as the file was opened; on a page just loaded, the Files list can
come before the page's first app state, so the file was kept as of no study,
and once the app state named the study, the Agent was never told the file
open (`previewing()` stayed null; the browser leg of CI failed on it,
waiting a minute for `notes.txt`). A text file is now kept as of the study
the server read it from, as it says with the text: also right where an
address names another study, or another study is opened as the page loads.
No file is said while no study is open.
"""

from __future__ import annotations

import json
import urllib.request

import pytest

pytest.importorskip("playwright.sync_api")

# The page's answers named here held back until the test lets each go, as a
# slow machine or network has them.
HELD = """(() => {
    const real = window.fetch.bind(window);
    const held = {};
    const gate = (name) => held[name] || (held[name] = (() => {
        let release;
        const released = new Promise((done) => { release = done; });
        return {released, release};
    })());
    window.__release = (name) => gate(name).release();
    window.fetch = async (input, ...rest) => {
        const url = String(typeof input === 'string' ? input : input.url);
        for (const name of window.__held || []) {
            if (url.includes(name)) await gate(name).released;
        }
        return real(input, ...rest);
    };
})();"""


def _page(browser, held):
    context = browser.new_context(viewport={"width": 1440, "height": 900})
    context.add_init_script(f"window.__held = {json.dumps(held)};")
    context.add_init_script(HELD)
    page = context.new_page()
    page.set_default_timeout(60000)
    page.errors = []
    page.on("pageerror", lambda error: page.errors.append(str(error)))
    return page


def _reading(page, path):
    page.wait_for_function("() => window.FastMDXFrame && window.FastMDXFrame.previewPath")
    page.evaluate("(path) => window.FastMDXFrame.previewPath(path)", path)
    page.wait_for_function("(path) => document.getElementById('side-preview-name')"
                           ".textContent === path", arg=path)


def _shown_study(page):
    return page.evaluate("""() => {
        const dash = window.FastMDXDashboard && window.FastMDXDashboard.state;
        return (dash && dash.appState && dash.appState.active_run) || null;
    }""")


def test_a_file_opened_before_the_app_state_is_the_one_read(tmp_path) -> None:
    from playwright.sync_api import sync_playwright

    from fastmdxplora.gui.server import start_dashboard_session
    from tests.test_the_drawing_scripts_run_in_a_browser import _write_study

    root = _write_study(tmp_path / "study")
    (root / "notes.txt").write_text("The ligand is benzamidine.\n", encoding="utf-8")
    session = start_dashboard_session(output=str(root), host="127.0.0.1", port=0)
    try:
        with sync_playwright() as pw:
            browser = pw.chromium.launch(args=["--enable-unsafe-swiftshader"])
            page = _page(browser, ["/api/app-state"])
            page.goto(session.url + "#files", wait_until="domcontentloaded")
            _reading(page, "notes.txt")
            before = _shown_study(page)
            page.evaluate("window.__release('/api/app-state')")
            page.wait_for_function("() => window.FastMDXDashboard.state.appState.active_run")
            read = page.evaluate("window.FastMDXFrame.previewing()")
            errors = page.errors
            browser.close()
    finally:
        session.server.shutdown()
    # Opened before the page knew the study, as on a slow load.
    assert before is None
    assert read == "notes.txt"
    assert errors == []


def test_a_file_of_the_study_before_is_not_said_open_in_the_next(tmp_path) -> None:
    """The page opened on an address naming another study: the Files list,
    read first, was of the study the server had open, and its file must not
    be said as open in the study opened in its place."""
    from playwright.sync_api import sync_playwright

    from fastmdxplora.gui.server import start_dashboard_session
    from fastmdxplora.gui.workspace import study_key
    from tests.test_the_workspace_says_its_studies import _study

    workspace = tmp_path / "workspace"
    first = _study(workspace / "ubiquitin", system="1UBQ", started="2026-09-01T10:00:00+00:00")
    second = _study(workspace / "chignolin", system="5AWL", started="2026-09-02T10:00:00+00:00")
    (first / "notes.txt").write_text("The first study's notes.\n", encoding="utf-8")
    (second / "notes.txt").write_text("The second study's notes.\n", encoding="utf-8")
    session = start_dashboard_session(output=str(first), host="127.0.0.1", port=0)
    base = session.url.rstrip("/")

    def open_now():
        with urllib.request.urlopen(base + "/api/app-state", timeout=30) as answer:
            return json.loads(answer.read()).get("active_run")

    try:
        with sync_playwright() as pw:
            browser = pw.chromium.launch(args=["--enable-unsafe-swiftshader"])
            page = _page(browser, ["/api/explore/switch", "/api/app-state"])
            page.goto(f"{base}/?study={study_key(second)}#files", wait_until="domcontentloaded")
            _reading(page, "notes.txt")
            # The text is read after the name is shown: waited for.
            page.wait_for_function("() => /study's notes/.test("
                                   "document.getElementById('side-preview-body').textContent)")
            text = page.text_content("#side-preview-body")
            page.evaluate("window.__release('/api/explore/switch')")
            for _ in range(300):
                if open_now() == str(second.resolve()):
                    break
                page.wait_for_timeout(100)
            page.evaluate("window.__release('/api/app-state')")
            page.wait_for_function("(second) => window.FastMDXDashboard.state.appState"
                                   ".active_run === second", arg=str(second.resolve()))
            read = page.evaluate("window.FastMDXFrame.previewing()")
            errors = page.errors
            browser.close()
    finally:
        session.server.shutdown()
    assert "first study" in text
    assert read is None
    assert errors == []


def test_a_text_is_said_of_the_study_it_was_read_from(tmp_path, monkeypatch) -> None:
    """Another study opened while a file was being read: the text and the
    study said with it are of the one it was read from."""
    import threading
    from pathlib import Path

    from fastmdxplora.gui import agent_panel
    from fastmdxplora.gui.server import start_dashboard_session
    from tests.test_the_workspace_says_its_studies import _study

    workspace = tmp_path / "workspace"
    first = _study(workspace / "ubiquitin", system="1UBQ", started="2026-09-01T10:00:00+00:00")
    second = _study(workspace / "chignolin", system="5AWL", started="2026-09-02T10:00:00+00:00")
    (first / "notes.txt").write_text("The first study's notes.\n", encoding="utf-8")
    (second / "notes.txt").write_text("The second study's notes.\n", encoding="utf-8")
    reading, read = threading.Event(), threading.Event()
    real = agent_panel.read_text_file

    def slowly(target, **kwargs):
        reading.set()
        read.wait(30)
        return real(target, **kwargs)

    monkeypatch.setattr(agent_panel, "read_text_file", slowly)
    session = start_dashboard_session(output=str(first), host="127.0.0.1", port=0)
    base = session.url.rstrip("/")
    answer: dict = {}

    def ask():
        with urllib.request.urlopen(base + "/api/file-text?path=notes.txt", timeout=60) as reply:
            answer.update(json.loads(reply.read()))

    try:
        asking = threading.Thread(target=ask)
        asking.start()
        assert reading.wait(30)
        switch = urllib.request.Request(
            base + "/api/explore/switch", data=json.dumps({"folder": str(second)}).encode(),
            headers={"Content-Type": "application/json", "Origin": base}, method="POST")
        with urllib.request.urlopen(switch, timeout=60) as reply:
            assert json.loads(reply.read())["ok"]
        read.set()
        asking.join(60)
    finally:
        read.set()
        session.server.shutdown()
    assert "first study" in answer["text"]
    assert Path(answer["study"]).resolve() == first.resolve()


def test_a_file_opened_from_a_list_of_the_study_before_is_of_the_one_open(tmp_path) -> None:
    """Another study opened elsewhere, the page told, its list not yet read
    again: a file opened from the old list is read from the study open, and
    is said as open in it (a file kept as of the list's study, as one version
    of this fix had it, was never said)."""
    from playwright.sync_api import sync_playwright

    from fastmdxplora.gui.server import start_dashboard_session
    from tests.test_the_workspace_says_its_studies import _study

    workspace = tmp_path / "workspace"
    first = _study(workspace / "ubiquitin", system="1UBQ", started="2026-09-01T10:00:00+00:00")
    second = _study(workspace / "chignolin", system="5AWL", started="2026-09-02T10:00:00+00:00")
    (first / "notes.txt").write_text("The first study's notes.\n", encoding="utf-8")
    (second / "notes.txt").write_text("The second study's notes.\n", encoding="utf-8")
    session = start_dashboard_session(output=str(first), host="127.0.0.1", port=0)
    base = session.url.rstrip("/")
    try:
        with sync_playwright() as pw:
            browser = pw.chromium.launch(args=["--enable-unsafe-swiftshader"])
            page = _page(browser, [])
            page.goto(base + "/#files", wait_until="domcontentloaded")
            page.wait_for_function("() => window.FastMDXDashboard.state.appState.active_run")
            page.evaluate("() => window.FastMDXFrame.showTab('files')")
            page.wait_for_selector("#side-files .side-file[data-path='notes.txt']",
                                   state="attached")
            # The list read again is held: the one shown is the first study's.
            page.evaluate("() => { window.__held = ['/api/artifacts']; }")
            switch = urllib.request.Request(
                base + "/api/explore/switch", data=json.dumps({"folder": str(second)}).encode(),
                headers={"Content-Type": "application/json", "Origin": base}, method="POST")
            with urllib.request.urlopen(switch, timeout=60) as reply:
                assert json.loads(reply.read())["ok"]
            page.wait_for_function("(second) => window.FastMDXDashboard.state.appState"
                                   ".active_run === second", arg=str(second.resolve()))
            _reading(page, "notes.txt")
            page.wait_for_function("() => /second study/.test("
                                   "document.getElementById('side-preview-body').textContent)")
            read = page.evaluate("window.FastMDXFrame.previewing()")
            page.evaluate("window.__release('/api/artifacts')")
            errors = page.errors
            browser.close()
    finally:
        session.server.shutdown()
    assert read == "notes.txt"
    assert errors == []


def test_a_file_opened_while_no_study_is_open_is_not_said(tmp_path) -> None:
    """The study closed (its data stale), the list not yet read again: a file
    opened from the old list reads as not found, and was said as open, with
    no study, or once the study was open again."""
    from playwright.sync_api import sync_playwright

    from fastmdxplora.gui.server import start_dashboard_session
    from tests.test_the_workspace_says_its_studies import _study

    first = _study(tmp_path / "workspace" / "ubiquitin", system="1UBQ",
                   started="2026-09-01T10:00:00+00:00")
    (first / "notes.txt").write_text("The first study's notes.\n", encoding="utf-8")
    session = start_dashboard_session(output=str(first), host="127.0.0.1", port=0)
    try:
        with sync_playwright() as pw:
            browser = pw.chromium.launch(args=["--enable-unsafe-swiftshader"])
            page = _page(browser, [])
            page.goto(session.url.rstrip("/") + "/#files", wait_until="domcontentloaded")
            page.wait_for_function("() => window.FastMDXDashboard.state.appState.active_run")
            page.evaluate("() => window.FastMDXFrame.showTab('files')")
            page.wait_for_selector("#side-files .side-file[data-path='notes.txt']",
                                   state="attached")
            page.evaluate("() => { window.__held = ['/api/artifacts']; }")
            session.runtime.data_stale = True
            page.wait_for_function("() => !window.FastMDXDashboard.state.appState.active_run")
            _reading(page, "notes.txt")
            # Its answer come (not found), before the study is open again.
            page.wait_for_selector("#side-preview-body .side-preview-note")
            closed = page.evaluate("window.FastMDXFrame.previewing()")
            session.runtime.data_stale = False
            page.wait_for_function("() => window.FastMDXDashboard.state.appState.active_run")
            back = page.evaluate("window.FastMDXFrame.previewing()")
            page.evaluate("window.__release('/api/artifacts')")
            errors = page.errors
            browser.close()
    finally:
        session.server.shutdown()
    assert (closed, back) == (None, None)
    assert errors == []


def test_a_file_read_from_the_next_study_is_not_said_in_the_one_before(tmp_path) -> None:
    """The page not yet told its study, another opened as a file of the list
    shown was opened: the reader shows that other study's file, said as open
    in it, and not in the first study once that is open again (a version of
    this fix, keeping the list's study, said it there)."""
    from playwright.sync_api import sync_playwright

    from fastmdxplora.gui.server import start_dashboard_session
    from tests.test_the_workspace_says_its_studies import _study

    workspace = tmp_path / "workspace"
    first = _study(workspace / "ubiquitin", system="1UBQ", started="2026-09-01T10:00:00+00:00")
    second = _study(workspace / "chignolin", system="5AWL", started="2026-09-02T10:00:00+00:00")
    (first / "notes.txt").write_text("The first study's notes.\n", encoding="utf-8")
    (second / "notes.txt").write_text("The second study's notes.\n", encoding="utf-8")
    session = start_dashboard_session(output=str(first), host="127.0.0.1", port=0)
    base = session.url.rstrip("/")

    def switch(folder):
        asked = urllib.request.Request(
            base + "/api/explore/switch", data=json.dumps({"folder": str(folder)}).encode(),
            headers={"Content-Type": "application/json", "Origin": base}, method="POST")
        with urllib.request.urlopen(asked, timeout=60) as reply:
            assert json.loads(reply.read())["ok"]

    try:
        with sync_playwright() as pw:
            browser = pw.chromium.launch(args=["--enable-unsafe-swiftshader"])
            page = _page(browser, ["/api/app-state"])
            page.goto(base + "/#files", wait_until="domcontentloaded")
            page.wait_for_function("() => window.FastMDXFrame && window.FastMDXFrame.showTab")
            page.evaluate("() => window.FastMDXFrame.showTab('files')")
            page.wait_for_selector("#side-files .side-file[data-path='notes.txt']",
                                   state="attached")
            switch(second)
            _reading(page, "notes.txt")
            page.wait_for_function("() => /second study/.test("
                                   "document.getElementById('side-preview-body').textContent)")
            page.evaluate("window.__release('/api/app-state')")
            page.wait_for_function("(second) => window.FastMDXDashboard.state.appState"
                                   ".active_run === second", arg=str(second.resolve()))
            in_second = page.evaluate("window.FastMDXFrame.previewing()")
            switch(first)
            page.wait_for_function("(first) => window.FastMDXDashboard.state.appState"
                                   ".active_run === first", arg=str(first.resolve()))
            in_first = page.evaluate("window.FastMDXFrame.previewing()")
            errors = page.errors
            browser.close()
    finally:
        session.server.shutdown()
    # The second study's file, said while it is open, and not in the first.
    assert (in_second, in_first) == ("notes.txt", None)
    assert errors == []


LATE_TEXT = """(() => {
    const real = window.fetch.bind(window);
    let asked = 0, release;
    const released = new Promise((done) => { release = done; });
    window.__releaseFirstText = () => release();
    window.fetch = async (input, ...rest) => {
        const url = String(typeof input === 'string' ? input : input.url);
        if (url.includes('/api/file-text') && ++asked === 1) {
            // Read by the server at once; its answer held.
            const answer = await real(input, ...rest);
            await released;
            return answer;
        }
        return real(input, ...rest);
    };
})();"""


def test_a_text_answered_late_shows_and_says_nothing(tmp_path) -> None:
    """A file opened, another study opened elsewhere, the file opened again:
    the first answer, late, drew the first study's text over the second's,
    under a file said as open in the second."""
    from playwright.sync_api import sync_playwright

    from fastmdxplora.gui.server import start_dashboard_session
    from tests.test_the_workspace_says_its_studies import _study

    workspace = tmp_path / "workspace"
    first = _study(workspace / "ubiquitin", system="1UBQ", started="2026-09-01T10:00:00+00:00")
    second = _study(workspace / "chignolin", system="5AWL", started="2026-09-02T10:00:00+00:00")
    (first / "notes.txt").write_text("The first study's notes.\n", encoding="utf-8")
    (second / "notes.txt").write_text("The second study's notes.\n", encoding="utf-8")
    session = start_dashboard_session(output=str(first), host="127.0.0.1", port=0)
    base = session.url.rstrip("/")
    try:
        with sync_playwright() as pw:
            browser = pw.chromium.launch(args=["--enable-unsafe-swiftshader"])
            page = _page(browser, [])
            page.context.add_init_script(LATE_TEXT)
            page.goto(base + "/#files", wait_until="domcontentloaded")
            page.wait_for_function("() => window.FastMDXDashboard.state.appState.active_run")
            page.evaluate("() => window.FastMDXFrame.showTab('files')")
            page.wait_for_selector("#side-files .side-file[data-path='notes.txt']",
                                   state="attached")
            # The page told nothing more until the test says so.
            page.evaluate("() => { window.__held = ['/api/app-state']; }")
            page.evaluate("() => window.FastMDXFrame.previewPath('notes.txt')")
            page.wait_for_function("() => window.__releaseFirstText")
            switch = urllib.request.Request(
                base + "/api/explore/switch", data=json.dumps({"folder": str(second)}).encode(),
                headers={"Content-Type": "application/json", "Origin": base}, method="POST")
            page.wait_for_timeout(1000)
            with urllib.request.urlopen(switch, timeout=60) as reply:
                assert json.loads(reply.read())["ok"]
            page.evaluate("() => window.FastMDXFrame.previewPath('notes.txt')")
            page.wait_for_function("() => /second study/.test("
                                   "document.getElementById('side-preview-body').textContent)")
            page.evaluate("() => window.__releaseFirstText()")
            page.wait_for_timeout(1000)
            shown = page.text_content("#side-preview-body")
            page.evaluate("window.__release('/api/app-state')")
            page.wait_for_function("(second) => window.FastMDXDashboard.state.appState"
                                   ".active_run === second", arg=str(second.resolve()))
            said = page.evaluate("window.FastMDXFrame.previewing()")
            errors = page.errors
            browser.close()
    finally:
        session.server.shutdown()
    assert "second study" in shown and "first study" not in shown
    assert said == "notes.txt"
    assert errors == []
