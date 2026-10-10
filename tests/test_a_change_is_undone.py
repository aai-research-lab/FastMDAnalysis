"""A change that can be undone is offered Undo in the notice after it.

Tags and a note kept on All studies, and a view saved in the Viewer or
forgotten there, were final: a slip took the person retyping what was
there, or setting the view up again by hand. The notice that says the
change now carries Undo for a while, which puts back what was there before;
a view saved over another is that one again, a new one goes.
"""

from __future__ import annotations

import json
import urllib.request

import pytest

from tests import viewer_hooks as hooks
from tests.test_the_workspace_says_its_studies import _study

pytest.importorskip("playwright.sync_api")

CAMERA = {"position": [1.0, 2.0, 30.0], "target": [1.0, 2.0, 3.0], "up": [0.0, 1.0, 0.0],
          "radius": 12.5, "fov": 0.785, "mode": "perspective"}


@pytest.fixture(scope="module")
def session(tmp_path_factory):
    from fastmdxplora.gui.server import start_dashboard_session
    from fastmdxplora.study_tags import set_tags

    workspace = tmp_path_factory.mktemp("undo")
    first = _study(workspace / "wild_type", started="2026-09-02T10:00:00+00:00")
    _study(workspace / "mutant", system="1UBQ", started="2026-09-03T10:00:00+00:00")
    set_tags(first, ["reference"], "the wild type, 300 K")
    started = start_dashboard_session(output=str(first), host="127.0.0.1", port=0)
    started.workspace, started.first = workspace, first
    yield started
    started.server.shutdown()


@pytest.fixture(scope="module")
def browser():
    from playwright.sync_api import sync_playwright

    with sync_playwright() as pw:
        launched = pw.chromium.launch(args=["--enable-unsafe-swiftshader"])
        yield launched
        launched.close()


def _page(browser, session, where):
    page = browser.new_context(viewport={"width": 1400, "height": 1000}).new_page()
    page.set_default_timeout(60000)
    page.errors = []
    page.on("pageerror", lambda error: page.errors.append(str(error)))
    page.goto(session.url + where, wait_until="domcontentloaded")
    return page


def _views(session):
    with urllib.request.urlopen(session.url.rstrip("/") + "/api/views", timeout=30) as answer:
        return json.loads(answer.read())["views"]


def _save(session, name, view):
    base = session.url.rstrip("/")
    request = urllib.request.Request(base + "/api/views", data=json.dumps(
        {"action": "save", "name": name, "view": view}).encode(),
        headers={"Content-Type": "application/json", "Origin": base}, method="POST")
    with urllib.request.urlopen(request, timeout=30) as answer:
        assert json.loads(answer.read())["ok"]


def _undo(page, said):
    page.wait_for_selector(f"#dashboard-toast.show:has-text('{said}') .dashboard-toast-action")
    page.click("#dashboard-toast .dashboard-toast-action")


def test_tags_and_a_note_kept_are_undone(browser, session) -> None:
    from fastmdxplora.study_tags import tags_of

    before = tags_of(session.first)
    page = _page(browser, session, "#studies")
    page.evaluate(f"() => window.FastMDXStudies.load({str(session.workspace)!r})")
    card = '.study-card[data-path$="wild_type"]'
    page.wait_for_selector(card)
    page.click(f"{card} .study-tag-edit")
    page.fill(f"{card} .study-tag-input", "slip")
    page.fill(f"{card} .study-note-input", "a note typed over the other")
    page.click(f"{card} .study-tag-editor button:has-text('Save')")
    page.wait_for_selector(f"{card} .study-tag:has-text('slip')")
    kept = tags_of(session.first)
    _undo(page, "Kept the tags and note")
    page.wait_for_selector("#dashboard-toast:has-text('back as they were')")
    shown = page.text_content(card)
    errors = page.errors
    page.context.close()
    assert "slip" in kept["tags"]
    assert tags_of(session.first) == before
    assert "slip" not in shown and "the wild type, 300 K" in shown
    assert errors == []


def test_a_view_forgotten_is_undone(browser, session) -> None:
    first = {"camera": CAMERA, "frame": 3}
    _save(session, "one", first)
    page = _page(browser, session, "#viewer")
    page.wait_for_selector("#viewer-views option[value='one']", state="attached")
    hooks.tool(page, "side-saved")
    page.select_option("#viewer-views", "one")
    page.click("#viewer-view-forget")
    page.wait_for_function("() => !document.querySelector(\"#viewer-views option[value='one']\")")
    gone = _views(session)
    _undo(page, "Forgot the view one")
    page.wait_for_selector("#viewer-views option[value='one']", state="attached")
    back = _views(session)
    errors = page.errors
    page.context.close()
    assert gone == []
    assert [(view["name"], view["frame"], view["camera"]["radius"]) for view in back] == [("one", 3, 12.5)]
    assert errors == []
