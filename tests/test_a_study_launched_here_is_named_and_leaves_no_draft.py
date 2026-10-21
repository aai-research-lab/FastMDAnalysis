"""A study a run started here: named by its folder until it has a manifest,
not by the study shown before it; and a study started from the Config
Builder leaves no draft of itself for the page opened next (CI-774)."""

from __future__ import annotations

import json
import urllib.request

from tests.test_open_studies_have_tabs import _builder, _open, browser, site  # noqa: F401


def _switch_elsewhere(session, folder):
    """Another window (or an AI app) opens ``folder``."""
    base = session.url.split("#")[0].rstrip("/")
    asked = urllib.request.Request(base + "/api/explore/switch", method="POST",
                                   data=json.dumps({"folder": str(folder)}).encode(),
                                   headers={"content-type": "application/json"})
    with urllib.request.urlopen(asked, timeout=30) as answer:
        assert json.loads(answer.read()).get("ok") is not False


def _shows(page, key, timeout=30000):
    page.wait_for_function("(key) => window.FastMDXDashboard.state.studyShown === key",
                           arg=key, timeout=timeout)


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


def test_a_study_with_no_manifest_is_not_named_by_the_study_before(browser, site):  # noqa: F811
    """A folder with no manifest (a run just started) opened after a study
    named by its system: it is named by its folder, as its tab and its
    status name it, not by the other study's system (L2: "BETA was opened
    elsewhere", said of a run's study)."""
    from fastmdxplora.gui.workspace import study_key

    bare = site["root"] / "bare"
    bare.mkdir()
    (bare / "exploration.yml").write_text("systems: []\n", encoding="utf-8")
    page = _open(browser, site)
    try:
        _switch_elsewhere(site["session"], site["studies"][1])
        _shows(page, site["keys"][1])
        page.wait_for_function("() => document.getElementById('topbar-run-title').textContent === 'BETA'")
        _switch_elsewhere(site["session"], bare)
        _shows(page, study_key(bare))
        page.wait_for_timeout(1500)
        title = page.text_content("#topbar-run-title")
        system = page.text_content("#topbar-run-id")
        tab = page.text_content(f'.study-tab-link[data-study-key="{study_key(bare)}"]')
        errors = list(page.errors)
    finally:
        page.context.close()
    assert title == "bare" and "bare" in tab, (title, tab)
    assert "BETA" not in system, system
    assert errors == []
