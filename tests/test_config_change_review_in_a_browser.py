"""Choose Config changes through the real workspace and keep the result."""
from copy import deepcopy
import json
import urllib.request

import pytest
import yaml

sync_api = pytest.importorskip("playwright.sync_api")
expect = sync_api.expect

BASE = {"systems": [{"id": "ubq", "system": "1UBQ"}], "include_phase": ["simulation"],
        "simulation": {"temperature_K": 300, "duration_ns": 2}}
LIVE = '#agent-thread .agent-study:not(.is-replaced):not([hidden])'


def _entry(config, version):
    from fastmdxplora.gui.plan import plan_of, sourced

    return {"role": "agent", "kind": "config", "config": config,
            "yaml": yaml.safe_dump(config), "version": version,
            "plan": sourced(plan_of(config), config), "attempts": [], "cycles": 1}


def _request(session, path, body=None):
    request = urllib.request.Request(session.url + path,
        data=json.dumps(body).encode() if body is not None else None,
        headers={"Origin": session.url, "Content-Type": "application/json"})
    with urllib.request.urlopen(request, timeout=10) as response:
        return json.load(response)


@pytest.fixture(scope="module")
def browser():
    from playwright.sync_api import sync_playwright

    with sync_playwright() as pw:
        browser = pw.chromium.launch(args=["--enable-unsafe-swiftshader"])
        yield browser
        browser.close()


@pytest.fixture
def site(tmp_path, monkeypatch):
    from fastmdxplora.gui.server import start_dashboard_session

    monkeypatch.setenv("FASTMDXPLORA_CONFIG_DIR", str(tmp_path / "config"))
    session = start_dashboard_session(output=str(tmp_path / "workspace"), host="127.0.0.1", port=0)
    yield session
    session.stop()


def _open(browser, site, before=None, after=None, width=1400, theme="dark"):
    from fastmdxplora.gui.plan import changes_between

    before = deepcopy(before or BASE)
    after = deepcopy(after or before)
    if after == before:
        after["simulation"].update(temperature_K=330, duration_ns=10)
    second = _entry(after, 2)
    second["changes"] = changes_between(before, after)
    _request(site, "/api/agent/conversation", {"entries": [
        {"role": "user", "text": "A study of ubiquitin."}, _entry(before, 1),
        {"role": "user", "text": "Change the conditions."}, second]})
    context = browser.new_context(viewport={"width": width, "height": 900})
    context.add_init_script(f"localStorage.setItem('fmx.theme', '{theme}');")
    page = context.new_page()
    page.set_default_timeout(7000)
    errors = []
    page.on("pageerror", lambda error: errors.append(str(error)))
    page.errors = errors
    page.goto(site.url + "#agent", wait_until="domcontentloaded")
    page.locator(LIVE).get_by_role("button", name="Review each change", exact=True).click()
    page.locator('[data-review-setting="simulation.temperature_K"]').wait_for()
    return page


def _choose(page, name, choice):
    page.locator(f'{LIVE} [data-review-setting="{name}"] [data-choice="{choice}"]').click()


def _apply(page):
    page.locator(LIVE).get_by_role("button", name="Use reviewed version", exact=True).click()
    page.wait_for_function("() => document.querySelectorAll('#agent-thread .agent-study:not([hidden])').length === 3")
    page.wait_for_function("""async () => {
        const data = await (await fetch('/api/agent/conversation')).json();
        return (data.entries || []).filter(e => e.kind === 'config').length === 3;
    }""")


def test_mixed_choices_become_a_new_version_and_survive_reload(browser, site):
    page = _open(browser, site)
    try:
        assert page.locator(LIVE).get_by_role("button", name="Use reviewed version").is_disabled()
        _choose(page, "simulation.temperature_K", "reject")
        _choose(page, "simulation.duration_ns", "accept")
        _apply(page)
        configs = [e for e in _request(site, "/api/agent/conversation")["entries"] if e.get("kind") == "config"]
        assert [e["config"]["simulation"] for e in configs] == [
            {"temperature_K": 300, "duration_ns": 2}, {"temperature_K": 330, "duration_ns": 10},
            {"temperature_K": 300, "duration_ns": 10}]
        assert configs[-1]["review"] == {"base": 1, "proposal": 2,
            "accepted": ["simulation.duration_ns"], "rejected": ["simulation.temperature_K"]}
        page.reload(wait_until="domcontentloaded")
        page.locator(LIVE + ' [data-role="version"]').filter(has_text="3").wait_for()
        assert page.locator('#agent-thread .agent-study.is-replaced [data-role="run"]:visible').count() == 0
        page.locator(LIVE + ' [data-role="load"]').click()
        page.locator('#run-settings [data-setting="temperature_K"] input').wait_for()
        assert page.input_value('#run-settings [data-setting="temperature_K"] input') == "300"
        assert page.input_value('#run-settings [data-setting="duration_ns"] input') == "10"
        assert page.errors == []
    finally:
        page.context.close()


def test_run_and_next_message_use_the_reviewed_config(browser, site):
    page = _open(browser, site)
    try:
        _choose(page, "simulation.temperature_K", "accept")
        _choose(page, "simulation.duration_ns", "reject")
        _apply(page)
        runs, prompts = [], []
        def run(route):
            runs.append(route.request.post_data_json)
            route.fulfill(status=200, content_type="application/json",
                          body=json.dumps({"ok": False, "error": "Run captured for this check."}))
        def prompt(route):
            prompts.append(route.request.post_data_json)
            route.fulfill(status=200, content_type="application/json",
                          body=json.dumps({"ok": False, "answer": "The reviewed conditions.", "cites": []}))
        page.route("**/api/agent/run", run)
        page.route("**/api/agent/propose*", prompt)
        page.locator(LIVE + ' [data-role="run"]').click()
        page.get_by_text("Run captured for this check.", exact=True).first.wait_for()
        assert runs[0]["config"]["simulation"] == {"temperature_K": 330, "duration_ns": 2}
        page.fill("#agent-request", "Describe these conditions.")
        page.keyboard.press("Enter")
        page.get_by_text("The reviewed conditions.", exact=True).wait_for()
        assert yaml.safe_load(prompts[0]["current_config"])["simulation"] == {"temperature_K": 330, "duration_ns": 2}
        assert page.errors == []
    finally:
        page.context.close()


def test_invalid_combination_creates_no_version_and_can_be_corrected(browser, site):
    before = deepcopy(BASE)
    before["analysis"] = {"include": ["rmsd"]}
    after = deepcopy(BASE)
    after["analysis"] = {"exclude": ["rg"]}
    after["simulation"]["temperature_K"] = 330
    page = _open(browser, site, before, after)
    try:
        page.locator(LIVE).get_by_role("button", name="Accept all", exact=True).click()
        _choose(page, "analysis.include", "reject")
        page.locator(LIVE).get_by_role("button", name="Use reviewed version").click()
        page.locator(LIVE + ' .agent-review-note').filter(has_text="mutually exclusive").wait_for()
        assert page.locator('#agent-thread .agent-study:not([hidden])').count() == 2
        assert page.locator(LIVE + ' [data-role="run"]').is_hidden()
        _choose(page, "analysis.include", "accept")
        _apply(page)
        config = [e for e in _request(site, "/api/agent/conversation")["entries"] if e.get("kind") == "config"][-1]["config"]
        assert config["analysis"] == {"exclude": ["rg"]}
        assert page.errors == []
    finally:
        page.context.close()


def test_pending_review_holds_the_message_and_cancel_restores_the_proposal(browser, site):
    page = _open(browser, site)
    try:
        requests = []
        page.route("**/api/agent/propose*", lambda route: (requests.append(route.request), route.abort()))
        _choose(page, "simulation.temperature_K", "reject")
        page.fill("#agent-request", "Run it.")
        page.keyboard.press("Enter")
        page.locator(LIVE + ' .agent-review-note').filter(has_text="before sending").wait_for()
        assert page.input_value("#agent-request") == "Run it." and requests == []
        page.locator(LIVE).get_by_role("button", name="Cancel review", exact=True).click()
        assert page.locator(LIVE + ' [data-role="run"]').is_visible()
        assert page.locator('#agent-thread .agent-study:not([hidden])').count() == 2
        expect(page.locator(LIVE).get_by_role("button", name="Review each change")).to_be_focused()
        assert page.errors == []
    finally:
        page.context.close()


def test_a_late_review_reply_cannot_enter_a_new_conversation(browser, site):
    page = _open(browser, site)
    try:
        page.locator(LIVE).get_by_role("button", name="Reject all", exact=True).click()
        held = []
        page.route("**/api/agent/review-config", lambda route: held.append(route))
        page.locator(LIVE).get_by_role("button", name="Use reviewed version").click()
        page.wait_for_function("() => document.querySelector('.agent-review-note').textContent.includes('Checking')")
        assert held
        page.click("#agent-new")
        page.wait_for_function("() => document.querySelectorAll('#agent-thread .agent-study:not([hidden])').length === 0")
        from fastmdxplora.gui.config_review import review_config
        checked = review_config(held[0].request.post_data_json)
        held[0].fulfill(status=200, content_type="application/json", body=json.dumps(checked))
        # A fetch after the held one provides an event-loop barrier, not a timed sleep.
        page.evaluate("async () => { await fetch('/api/agent/conversation'); await new Promise(requestAnimationFrame); }")
        assert page.locator('#agent-thread .agent-study:not([hidden])').count() == 0
        assert "Your choices" not in page.text_content("#agent-thread")
        assert page.errors == []
    finally:
        page.context.close()


@pytest.mark.parametrize("width", [390, 980, 1600])
@pytest.mark.parametrize("theme", ["light", "dark"])
def test_review_controls_fit_and_work_by_keyboard(browser, site, width, theme):
    page = _open(browser, site, width=width, theme=theme)
    try:
        focused = page.locator(LIVE + ' [data-review-setting="simulation.temperature_K"] [data-choice="accept"]')
        expect(focused).to_be_focused()
        page.keyboard.press("Space")
        assert focused.get_attribute("aria-pressed") == "true"
        page.keyboard.press("Tab")
        page.keyboard.press("Space")
        assert page.locator(LIVE + ' [data-review-setting="simulation.temperature_K"] [data-choice="reject"]').get_attribute("aria-pressed") == "true"
        bounds = page.locator(LIVE + ' [data-role="review"]').evaluate("""host => {
            const box = host.getBoundingClientRect();
            return [...host.querySelectorAll('fieldset, button')].filter(el => {
                const r = el.getBoundingClientRect();
                return r.left < box.left - 1 || r.right > box.right + 1;
            }).map(el => el.outerHTML);
        }""")
        assert bounds == []
        assert page.evaluate("document.documentElement.scrollWidth <= innerWidth")
        assert page.errors == []
    finally:
        page.context.close()


def test_alias_names_share_one_visible_choice_and_one_result(browser, site):
    before = deepcopy(BASE)
    before["simulation"].update(prepared_from="runs/old", setup_from="runs/old")
    after = deepcopy(before)
    after["simulation"].update(prepared_from="runs/new", setup_from="runs/new", temperature_K=330)
    after["decisions"] = {"simulation.setup_from": {"why": "Use the new preparation.", "source": "person"}}
    page = _open(browser, site, before, after)
    try:
        rows = page.locator(LIVE + ' [data-review-setting="simulation.prepared_from"]')
        assert rows.count() == 1
        assert page.locator(LIVE + ' [data-review-setting="simulation.setup_from"]').count() == 0
        assert "Equivalent setting names" in rows.text_content()
        assert "Use the new preparation." in rows.text_content()
        _choose(page, "simulation.temperature_K", "reject")
        _choose(page, "simulation.prepared_from", "accept")
        _apply(page)
        config = [e for e in _request(site, "/api/agent/conversation")["entries"] if e.get("kind") == "config"][-1]["config"]
        assert config["simulation"]["prepared_from"] == config["simulation"]["setup_from"] == "runs/new"
        assert config["simulation"]["temperature_K"] == 300
        assert config["decisions"]["simulation.setup_from"]["why"] == "Use the new preparation."
        assert page.errors == []
    finally:
        page.context.close()


def test_changed_alternatives_are_visible_before_choosing(browser, site):
    before = deepcopy(BASE)
    before["decisions"] = {"simulation.temperature_K": {"why": "Same reason.", "source": "person", "alternatives": [310]}}
    after = deepcopy(before)
    after["decisions"]["simulation.temperature_K"]["alternatives"] = [330]
    page = _open(browser, site, before, after)
    try:
        reason = page.locator(LIVE + ' [data-review-setting="simulation.temperature_K"] .agent-review-reasons').text_content()
        assert "alternatives: [310]" in reason and "alternatives: [330]" in reason
        _choose(page, "simulation.temperature_K", "accept")
        _apply(page)
        config = [e for e in _request(site, "/api/agent/conversation")["entries"] if e.get("kind") == "config"][-1]["config"]
        assert config["decisions"]["simulation.temperature_K"]["alternatives"] == [330]
        assert page.errors == []
    finally:
        page.context.close()


@pytest.mark.parametrize("full", [False, True])
def test_download_keeps_reviewed_settings_even_in_an_inactive_phase(browser, site, full):
    after = deepcopy(BASE)
    after["simulation"]["temperature_K"] = 330
    after["report"] = {"title": "An accepted title for later."}
    after["decisions"] = {"report.title": {"why": "The chosen title.", "source": "person"}}
    page = _open(browser, site, after=after)
    try:
        _choose(page, "simulation.temperature_K", "reject")
        _choose(page, "report.title", "accept")
        _apply(page)
        page.evaluate("window.showSaveFilePicker = undefined")
        page.locator(LIVE + ' [data-role="more"] > summary').click()
        if full:
            page.locator(LIVE + ' [data-role="full"]').check()
        with page.expect_download() as saved:
            page.locator(LIVE + ' [data-role="download"]').click()
        with open(saved.value.path(), encoding="utf-8") as stream:
            config = yaml.safe_load(stream)
        assert config["report"]["title"] == "An accepted title for later."
        assert config["decisions"]["report.title"]["why"] == "The chosen title."
        assert config["simulation"]["temperature_K"] == 300
        assert config["include_phase"] == ["simulation"]
        assert page.errors == []
    finally:
        page.context.close()
