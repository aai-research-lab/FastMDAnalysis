"""The Problems panel reads records, opens evidence and follows the study."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

pytest.importorskip("playwright.sync_api")


def _write_study(root: Path, why: str = "Only four independent samples.") -> Path:
    folder = root / "analysis" / "rmsd"
    folder.mkdir(parents=True)
    (root / "manifest.json").write_text(json.dumps({"phases": [
        {"name": "analysis", "status": "completed"}]}), encoding="utf-8")
    (root / "analysis" / "analysis_manifest.json").write_text(json.dumps({
        "results": {"rmsd": {"status": "success"}}}), encoding="utf-8")
    (folder / "options.json").write_text(json.dumps({"findings": {"mean": {
        "mean": 0.12, "unit": "nm", "effective_samples": 4,
        "not_a_measurement": why}}}), encoding="utf-8")
    return root


@pytest.fixture(scope="module")
def sessions(tmp_path_factory):
    from fastmdxplora.gui.server import start_dashboard_session

    base = tmp_path_factory.mktemp("problems")
    studies = {"first": _write_study(base / "first"),
               "second": _write_study(base / "second", "Second study's reason."),
               "empty": base / "empty"}
    studies["empty"].mkdir()
    sessions = {name: start_dashboard_session(output=str(path), host="127.0.0.1", port=0)
                for name, path in studies.items()}
    yield sessions, studies
    for session in sessions.values():
        session.stop()


@pytest.fixture(scope="module")
def browser():
    from playwright.sync_api import sync_playwright

    with sync_playwright() as pw:
        launched = pw.chromium.launch(args=["--enable-unsafe-swiftshader"])
        yield launched
        launched.close()


def _open(browser, session, *, width=1440, theme="dark"):
    context = browser.new_context(viewport={"width": width, "height": 900})
    context.add_init_script(
        f"localStorage.setItem('fmx.theme', '{theme}');"
        "localStorage.setItem('fmx.panelCollapsed', '0');")
    page = context.new_page()
    page.set_default_timeout(15000)
    page.goto(session.url + "#overview", wait_until="domcontentloaded")
    page.wait_for_function("() => !!window.FastMDXDashboard?.state.appState.mode")
    page.wait_for_function("() => document.getElementById('loading-screen').getAttribute('aria-hidden') === 'true'")
    if width <= 800:
        page.get_by_role("button", name="Study problems", exact=True).click()
    elif page.locator("#side-expand").is_visible():
        page.locator("#side-expand").click()
    page.locator("#problems-tab").click()
    return page


@pytest.mark.parametrize("width", [390, 900, 1440])
@pytest.mark.parametrize("theme", ["light", "dark"])
def test_the_recorded_reason_and_remedy_fit_the_panel(browser, sessions, width, theme):
    from playwright.sync_api import expect

    page = _open(browser, sessions[0]["first"], width=width, theme=theme)
    try:
        expect(page.locator("#problems-summary")).to_have_text("1 recorded problem")
        expect(page.locator("#problems-body")).to_contain_text("Only four independent samples.")
        expect(page.locator("#problems-body h4")).to_have_text("What would fix it")
        expect(page.locator("#problems-tab")).to_have_attribute("aria-selected", "true")
        assert page.evaluate("""() => {
            const pane = document.getElementById('problems-pane');
            const panel = pane.getBoundingClientRect();
            return panel.left >= -1 && panel.right <= innerWidth + 1
              && pane.scrollWidth <= pane.clientWidth + 1
              && [...pane.querySelectorAll('h3,h4,p,code,button')].every(el => {
                  const r = el.getBoundingClientRect();
                  return r.left >= panel.left - 1 && r.right <= panel.right + 1;
              });
        }""")
        assert page.locator("html").get_attribute("data-theme") == theme
    finally:
        page.context.close()


def test_evidence_opens_in_the_existing_files_reader(browser, sessions):
    from playwright.sync_api import expect

    page = _open(browser, sessions[0]["first"])
    try:
        page.get_by_role("button", name="Open analysis/rmsd/options.json in Files", exact=True).click()
        expect(page.locator('[data-side-tab="files"]')).to_have_attribute("aria-selected", "true")
        expect(page.locator("#side-preview-name")).to_contain_text("options.json")
        expect(page.locator("#side-preview-body")).to_contain_text("Only four independent samples.")
    finally:
        page.context.close()


def test_phone_evidence_and_escape_keep_the_panel_usable(browser, sessions):
    from playwright.sync_api import expect

    page = _open(browser, sessions[0]["first"], width=390)
    try:
        page.get_by_role("button", name="Open analysis/rmsd/options.json in Files", exact=True).click()
        expect(page.locator("#side-preview-body")).to_contain_text("Only four independent samples.")
        page.keyboard.press("Escape")
        expect(page.locator("#problems-mobile")).not_to_have_attribute("open", "")
        expect(page.get_by_role("button", name="Study problems", exact=True)).to_be_focused()
        page.get_by_role("button", name="Study problems", exact=True).click()
        expect(page.locator("#problems-summary")).to_have_text("1 recorded problem")
        page.set_viewport_size({"width": 1440, "height": 900})
        page.wait_for_function("() => !document.getElementById('problems-mobile').open")
        expect(page.locator(".app-shell > #side-panel")).to_have_count(1)
    finally:
        page.context.close()


def test_recorded_markup_is_shown_as_text(browser, sessions):
    from playwright.sync_api import expect

    page = _open(browser, sessions[0]["first"])
    try:
        page.route("**/api/problems?*", lambda route: route.fulfill(json={
            "ok": True, "study": page.evaluate("new URL(location.href).searchParams.get('study')"),
            "active_run": str(sessions[1]["first"]), "available": True, "notes": [],
            "items": [{"id": "x", "title": "<img src=x onerror=alert(1)>",
                       "reason": "<script>window.bad=true</script>", "fix": "<b>Read</b>"}]}))
        page.get_by_role("button", name="Refresh problems", exact=True).click()
        expect(page.locator("#problems-body")).to_contain_text("<script>window.bad=true</script>")
        assert page.locator("#problems-body img, #problems-body script, #problems-body b").count() == 0
        assert page.evaluate("window.bad === undefined")
    finally:
        page.context.close()


def test_a_failed_read_is_visible_and_can_be_retried(browser, sessions):
    from playwright.sync_api import expect

    page = _open(browser, sessions[0]["first"])
    try:
        expect(page.locator("#problems-summary")).to_have_text("1 recorded problem")
        page.route("**/api/problems?*", lambda route: route.fulfill(status=500, body="unavailable"))
        page.get_by_role("button", name="Refresh problems", exact=True).click()
        expect(page.locator("#problems-summary")).to_contain_text("Could not read problems")
        assert page.locator(".problem-item").count() == 0
        page.unroute("**/api/problems?*")
        page.get_by_role("button", name="Refresh problems", exact=True).click()
        expect(page.locator("#problems-summary")).to_have_text("1 recorded problem")
    finally:
        page.context.close()


def test_a_late_response_cannot_restore_the_previous_study(browser, sessions):
    from playwright.sync_api import expect

    page = _open(browser, sessions[0]["first"])
    try:
        expect(page.locator("#problems-summary")).to_have_text("1 recorded problem")
        # Let the HTTP request reach the old study, but hold its reply while
        # app-state moves on. The promise must remain harmless after that.
        page.evaluate("""() => {
            const original = window.fetch;
            window.__releaseProblems = null;
            window.fetch = (...args) => {
                const answer = original(...args);
                if (!String(args[0]).startsWith('/api/problems')) return answer;
                return answer.then(response => new Promise(resolve => {
                    window.__releaseProblems = () => resolve(response);
                }));
            };
        }""")
        page.get_by_role("button", name="Refresh problems", exact=True).click()
        page.wait_for_function("() => !!window.__releaseProblems")
        page.evaluate("""() => window.FastMDXDashboard.applyAppState({
            mode: 'run', active_run: '/different-study', status: 'idle'
        })""")
        assert page.locator(".problem-item").count() == 0
        page.evaluate("() => window.__releaseProblems()")
        expect(page.locator("#problems-summary")).to_contain_text("Reading the study's records")
        assert page.locator(".problem-item").count() == 0
        expect(page.get_by_role("button", name="Refresh problems", exact=True)).to_be_enabled()
    finally:
        page.context.close()


def test_no_study_does_not_claim_a_clean_bill_of_health(browser, sessions):
    from playwright.sync_api import expect

    page = _open(browser, sessions[0]["empty"])
    try:
        expect(page.locator("#problems-summary")).to_have_text("Open a study to read its problems.")
        assert page.locator(".problem-item").count() == 0
    finally:
        page.context.close()


def test_the_endpoint_refuses_to_read_a_different_study(sessions):
    from urllib.request import urlopen

    with urlopen(sessions[0]["first"].url + "/api/problems?study=wrong", timeout=30) as response:
        answer = json.load(response)
    assert answer["ok"] is False
    assert answer["reason"] == "The open study changed."
    assert "items" not in answer
