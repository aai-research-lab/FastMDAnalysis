"""Inspect actual audited atoms and accept the same retained plan through the GUI."""
from __future__ import annotations

import json

import pytest

from tests.test_setup_waits_for_acceptance import prepared as _prepared

prepared = _prepared

sync_api = pytest.importorskip("playwright.sync_api")
expect = sync_api.expect


@pytest.fixture(scope="module")
def browser():
    with sync_api.sync_playwright() as pw:
        instance = pw.chromium.launch(args=["--enable-unsafe-swiftshader"])
        yield instance
        instance.close()


@pytest.fixture
def site(prepared, monkeypatch, tmp_path):
    from fastmdxplora.gui.server import start_dashboard_session

    monkeypatch.setenv("FASTMDXPLORA_CONFIG_DIR", str(tmp_path / "config"))
    monkeypatch.setattr("fastmdxplora.gui.exploration.exploration_environment_error", lambda config: None)
    session = start_dashboard_session(output=str(prepared), exploration_root=tmp_path, host="127.0.0.1", port=0)
    calls = []
    monkeypatch.setattr(session.runtime, "_spawn", lambda command, root, url, **kwargs:
        (calls.append((command, root)), {"launched": True, "pid": 123})[1])
    yield {"session": session, "root": prepared, "calls": calls}
    session.stop()


def _open(browser, site, width=1440, theme="dark"):
    context = browser.new_context(viewport={"width": width, "height": 900})
    context.add_init_script(f"localStorage.setItem('fmx.theme', '{theme}')")
    page = context.new_page(); page.set_default_timeout(12000)
    page.errors = []; page.on("pageerror", lambda error: page.errors.append(str(error)))
    page.goto(site["session"].url + "#overview", wait_until="domcontentloaded")
    page.locator("#preparation-review-open").click()
    page.locator("#preparation-review-events button").first.wait_for()
    return page


@pytest.mark.parametrize("width,theme", [(390,"dark"),(390,"light"),(980,"dark"),(980,"light"),(1600,"dark"),(1600,"light")])
def test_observed_hydrogen_maps_to_the_actual_atom_in_both_snapshots(browser, site, width, theme):
    page = _open(browser, site, width, theme)
    try:
        expect(page.locator("#preparation-review-notice")).to_contain_text("does not certify")
        expect(page.locator("#preparation-review-choices")).to_contain_text("Recording a choice does not prove")
        page.get_by_role("button", name="add missing hydrogens", exact=True).click()
        page.locator("#preparation-change-canvas[data-ready='true']").wait_for()
        added = page.locator("#preparation-change-atoms button").filter(has_text="Added")
        expect(added).to_have_count(1)
        added.click()
        page.wait_for_function("window.FastMDXPreparationReview.engine?.scene.selected?.length === 1")
        atom = page.evaluate("""() => {const e=window.FastMDXPreparationReview.engine;return e.atoms(e.scene.selected)[0];}""")
        assert (atom["chain"], atom["resi"], atom["resn"], atom["atom"]) == ("A",1,"GLY","H")
        assert [atom["x"],atom["y"],atom["z"]] == pytest.approx([.7,2,3],abs=.005)
        page.locator("#preparation-before").click()
        page.wait_for_function("window.FastMDXPreparationReview.engine?.atomCount() === 4")
        page.locator("#preparation-after").click()
        page.wait_for_function("window.FastMDXPreparationReview.engine?.atomCount() === 5")
        assert page.evaluate("document.scrollingElement.scrollWidth <= innerWidth")
        assert page.locator(".preparation-review-box").bounding_box()["width"] <= width
        page.keyboard.press("Escape")
        expect(page.locator("#preparation-review-dialog")).to_be_hidden()
        expect(page.locator("#preparation-review-open")).to_be_focused()
        assert page.errors == []
    finally: page.context.close()


def test_acceptance_survives_reload_and_starts_only_the_existing_simulate_command(browser, site):
    page = _open(browser, site)
    try:
        expect(page.locator("#preparation-accept")).to_be_enabled()
        assert site["calls"] == []
        page.reload(wait_until="domcontentloaded")
        page.locator("#preparation-review-open").click()
        expect(page.locator("#preparation-accept")).to_be_enabled()
        assert site["calls"] == []
        page.locator("#preparation-accept").click()
        expect(page.locator("#preparation-review-status")).to_contain_text("Setup accepted")
        expect(page.locator("#preparation-accept")).to_be_disabled()
        assert len(site["calls"]) == 1
        command, root = site["calls"][0]
        assert command[3] == "simulate" and str(root) == str(site["root"])
        page.reload(wait_until="domcontentloaded")
        page.locator("#preparation-review-open").click()
        expect(page.locator("#preparation-accept")).to_be_disabled()
        assert len(site["calls"]) == 1 and page.errors == []
    finally: page.context.close()


def test_changed_preparation_is_refused_even_after_refresh(browser, site):
    page = _open(browser, site)
    try:
        expect(page.locator("#preparation-accept")).to_be_enabled()
        # Keep the browser's checked record stale deliberately. A live poll
        # may otherwise detect the change and correctly disable the button
        # before this test can exercise the server's acceptance refusal.
        observed = page.request.get(site["session"].url + "/api/preparation-review").json()
        page.route("**/api/preparation-review?*", lambda route:
            route.fulfill(content_type="application/json", body=json.dumps(observed)))
        target = site["root"] / "setup" / "system.xml"
        target.write_text("<System changed='true'/>")
        page.reload(wait_until="domcontentloaded")
        page.locator("#preparation-review-open").click()
        expect(page.locator("#preparation-accept")).to_be_enabled()
        page.locator("#preparation-accept").click()
        expect(page.locator("#preparation-review-status")).to_contain_text("changed")
        assert not site["calls"]
        page.unroute("**/api/preparation-review?*")
        page.locator("#preparation-review-refresh").click()
        expect(page.locator("#preparation-review-status")).to_contain_text("changed")
        expect(page.locator("#preparation-accept")).to_be_disabled()
        assert not site["calls"]
    finally: page.context.close()


def test_config_review_option_belongs_to_the_draft_and_is_cleared_on_load(browser, site):
    page = _open(browser, site)
    try:
        page.keyboard.press("Escape")
        page.evaluate("window.FastMDXDashboard.navigate('run')")
        page.wait_for_function("window.FastMDXRun?.ready")
        page.evaluate("window.FastMDXRun.setStart('structure')")
        page.locator("#run-setup-review").check()
        page.wait_for_function("window.FastMDXRun.state.reviewAfterSetup === true")
        page.reload(wait_until="domcontentloaded")
        page.wait_for_function("window.FastMDXRun?.ready")
        expect(page.locator("#run-setup-review")).to_be_checked()
        loaded = page.evaluate("""async () => (await (await fetch('/api/load-config',{method:'POST',
          headers:{'Content-Type':'application/json'},body:JSON.stringify({config:{systems:[{system:'1UBQ'}],
            include_phase:['setup','simulation']}})})).json()).state""")
        page.evaluate("state => window.FastMDXRun.applyLoadedState(state)", loaded)
        expect(page.locator("#run-setup-review")).not_to_be_checked()
        page.locator("#run-setup-review").check()
        page.get_by_role("button", name="Reset", exact=True).click()
        expect(page.locator("#run-setup-review")).not_to_be_checked()
        assert page.errors == []
    finally: page.context.close()


@pytest.mark.parametrize("switch_study", [False, True])
def test_late_observation_cannot_reopen_a_closed_or_previous_study(browser, site, tmp_path, switch_study):
    page = _open(browser, site)
    try:
        held = []
        query = page.evaluate("""async () => {
          const data = await (await fetch('/api/preparation-review')).json();
          return new URLSearchParams({study:data.study,version:data.version,event:'event-0002'}).toString();
        }""")
        answer = page.request.get(site["session"].url + "/api/preparation-change?" + query).json()
        page.route("**/api/preparation-change?*", lambda route: held.append(route))
        with page.expect_request("**/api/preparation-change?*"):
            page.get_by_role("button", name="add missing hydrogens", exact=True).click()
        if switch_study:
            other = tmp_path / "other"; other.mkdir()
            (other / "manifest.json").write_text('{"phases":[{"name":"setup","status":"ok"}]}')
            page.evaluate("""async folder => {
              const response=await fetch('/api/explore/switch',{method:'POST',headers:{'Content-Type':'application/json'},
                body:JSON.stringify({folder})});const answer=await response.json();
              if(!answer.ok)throw new Error(JSON.stringify(answer));window.FastMDXDashboard.studyOpened(answer,'overview');
            }""", str(other))
        else:
            page.keyboard.press("Escape")
        expect(page.locator("#preparation-review-dialog")).to_be_hidden()
        held[0].fulfill(content_type="application/json", body=json.dumps(answer))
        page.wait_for_function("window.FastMDXPreparationReview !== undefined")
        page.evaluate("() => new Promise(requestAnimationFrame)")
        assert page.evaluate("window.FastMDXPreparationReview.engine === null")
        assert page.locator("#preparation-change-canvas canvas").count() == 0
        assert page.errors == []
    finally: page.context.close()


def test_pending_acceptance_holds_workspace_switch_and_duplicate_clicks(browser, site):
    page = _open(browser, site)
    try:
        held = []
        page.route("**/api/accept-preparation", lambda route: held.append(route))
        with page.expect_request("**/api/accept-preparation"):
            page.locator("#preparation-accept").click()
        expect(page.locator("#preparation-accept")).to_be_disabled()
        page.locator("#preparation-accept").evaluate("button => button.click()")
        answer = page.evaluate("window.FastMDXDashboard.openStudyTab('another-study','overview')")
        assert not answer["ok"] and "acceptance" in answer["error"]
        assert len(held) == 1 and not site["calls"]
        held[0].fulfill(content_type="application/json", body=json.dumps({"ok":False,"error":"No GPU available."}))
        expect(page.locator("#preparation-review-status")).to_contain_text("No GPU available")
        expect(page.locator("#preparation-accept")).to_be_enabled()
        assert page.errors == []
    finally: page.context.close()


def test_preflight_failure_survives_actual_polling_until_refresh(browser, site, monkeypatch):
    monkeypatch.setattr("fastmdxplora.gui.exploration.exploration_environment_error",
                        lambda config: "Test backend is unavailable.")
    page = _open(browser, site)
    try:
        page.evaluate("""() => { window.reviewPolls = 0;
            window.FastMDXDashboard.on('app-state', () => { window.reviewPolls += 1; }); }""")
        page.locator("#preparation-accept").click()
        expect(page.locator("#preparation-review-status")).to_contain_text("Test backend is unavailable")
        polls = page.evaluate("window.reviewPolls")
        page.wait_for_function("n => window.reviewPolls >= n + 2", arg=polls)
        expect(page.locator("#preparation-review-status")).to_contain_text("Test backend is unavailable")
        assert not site["calls"]
        page.locator("#preparation-review-refresh").click()
        expect(page.locator("#preparation-review-status")).to_contain_text("waiting for your acceptance")
        assert not site["calls"] and page.errors == []
    finally: page.context.close()


def test_refusal_feedback_clears_when_the_prepared_files_change(browser, site, monkeypatch):
    monkeypatch.setattr("fastmdxplora.gui.exploration.exploration_environment_error",
                        lambda config: "Test backend is unavailable.")
    page = _open(browser, site)
    try:
        page.locator("#preparation-accept").click()
        expect(page.locator("#preparation-review-status")).to_contain_text("Test backend is unavailable")
        (site["root"] / "setup" / "system.xml").write_text("changed system")
        expect(page.locator("#preparation-review-status")).to_contain_text("changed")
        expect(page.locator("#preparation-accept")).to_be_disabled()
        assert not site["calls"] and page.errors == []
    finally: page.context.close()


def test_lost_acceptance_reply_requires_explicit_refresh_before_retry(browser, site):
    page = _open(browser, site)
    try:
        page.route("**/api/accept-preparation", lambda route: route.abort("failed"))
        page.locator("#preparation-accept").click()
        expect(page.locator("#preparation-review-status")).to_contain_text("server did not answer")
        expect(page.locator("#preparation-accept")).to_be_disabled()
        page.evaluate("window.dispatchEvent(new CustomEvent('dashboard:app-state', {detail:{}}))")
        expect(page.locator("#preparation-review-status")).to_contain_text("server did not answer")
        expect(page.locator("#preparation-accept")).to_be_disabled()
        page.locator("#preparation-review-refresh").click()
        expect(page.locator("#preparation-accept")).to_be_enabled()
        assert not site["calls"]
    finally: page.context.close()
