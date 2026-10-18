"""What is asked while a run the Agent started is being launched lands where
it was asked for (CI-774, second round).

Found by review K1: another study's conversation opened during a launch was
opened by the server at once, then the launch switched the GUI to the run's
folder, and the page opened whatever study the GUI had by then: the run's.
"""

from __future__ import annotations

import json

import pytest

from tests.test_a_conversation_follows_a_run_that_wrote_nothing import (  # noqa: F401
    _CONFIG,
    _browsers_closed,
    _gui,
    _page_on,
    _say,
    _slow_disk,
)

pytest.importorskip("playwright.sync_api")


def _other_study(workspace):
    other = workspace / "other"
    store = other / "agent" / "conversations"
    store.mkdir(parents=True)
    (other / "exploration.yml").write_text("systems: []\n", encoding="utf-8")
    (store / "conv-20261008-000000-000000.json").write_text(json.dumps(
        {"version": 3, "id": "conv-20261008-000000-000000",
         "entries": [{"eid": "o", "role": "user", "text": "about other"}]}), encoding="utf-8")
    return other


def test_another_study_opened_while_a_run_starts_is_the_one_opened(tmp_path, monkeypatch):
    from fastmdxplora.gui.workspace import study_key

    workspace, made, session = _gui(tmp_path, monkeypatch, taking=1.5)
    other = _other_study(workspace)
    asked = []
    try:
        page, errors, close = _page_on(session, [dict(_CONFIG)])
        _say(page, "chignolin")
        page.wait_for_selector("#agent-thread .agent-study:not([hidden]) [data-role=run]")
        page.wait_for_timeout(800)
        page.on("request", lambda r: asked.append(r.url.split("/api/")[-1])
                if "/api/agent/run" in r.url or "/api/agent/conversation/open" in r.url else None)
        _slow_disk(monkeypatch, 1.0)
        page.click("#agent-thread .agent-study:not([hidden]) [data-role=run]")
        page.click("#agent-conversations")
        page.wait_for_selector("#agent-conv-list:not([hidden])")
        page.evaluate("(where) => { window.__opened = window.FastMDXAgentPanel.open("
                      "'conv-20261008-000000-000000', where, false); }", str(other))
        said = page.evaluate("() => (document.querySelector("
                             "'#agent-conv-list .agent-conv-held') || {}).textContent || ''")
        page.wait_for_timeout(9000)
        status = page.evaluate("async () => (await (await fetch('/api/status')).json()).study")
        shown = page.evaluate("() => window.FastMDXDashboard.state.studyShown")
        address = page.evaluate("() => new URL(location.href).searchParams.get('study')")
        close()
    finally:
        session.server.shutdown()
    assert made
    assert status == shown == address == study_key(other), (status, shown, address)
    # The whole open waited for the launch: asked of the server only then.
    assert asked[:2] == ["agent/run", "agent/conversation/open"], asked
    assert said == "Opens once the run has started."
    assert errors == []
