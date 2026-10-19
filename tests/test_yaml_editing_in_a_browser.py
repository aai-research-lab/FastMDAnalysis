"""YAML edits remain a draft until their checked Builder conversion is applied."""
from __future__ import annotations

import json

import pytest
import yaml

from tests.test_the_builder_page_in_a_browser import site as builder_site, _open as builder_page

sync_api = pytest.importorskip("playwright.sync_api")
expect = sync_api.expect
site = builder_site


@pytest.fixture(autouse=True)
def no_scientific_start(monkeypatch):
    monkeypatch.setattr("fastmdxplora.gui.exploration.DashboardRuntime._spawn",
                        lambda *args, **kwargs: pytest.fail("YAML editing attempted to start a workflow."))


def document(site, temperature=310):
    return yaml.safe_dump({"systems":[{"system": str(site["root"] / "peptide.pdb")}],
        "include_phase":["setup","simulation"], "simulation":{"temperature_K":temperature}}, sort_keys=False)


def editor(site, width=1440, theme="dark"):
    context = site["browser"].new_context(viewport={"width":width,"height":900})
    context.add_init_script(f"localStorage.setItem('fmx.theme', '{theme}')")
    page = builder_page(site, context=context)
    page.wait_for_function("window.FastMDXRun.ready")
    page.set_default_timeout(12000)
    page.locator("#config-yaml-open").click()
    expect(page.locator("#config-yaml-dialog")).to_be_visible()
    return page


@pytest.mark.parametrize("width,theme", [(390,"dark"),(390,"light"),(980,"dark"),(980,"light"),(1600,"dark"),(1600,"light")])
def test_line_links_select_the_refused_yaml_and_work_in_both_themes(site, width, theme):
    page = editor(site, width, theme)
    try:
        text = document(site) + "simulation_typo: true\n"
        page.locator("#config-yaml-text").fill(text)
        page.locator("#config-yaml-validate").click()
        expect(page.locator("#config-yaml-status")).to_contain_text("Unknown top-level key")
        expected_line = len(text.splitlines())
        page.get_by_role("button", name="Line " + str(expected_line), exact=True).click()
        selected = page.evaluate("() => {const t=document.getElementById('config-yaml-text');return t.value.slice(t.selectionStart,t.selectionEnd)}")
        assert selected == "simulation_typo: true"
        expect(page.locator("#config-yaml-text")).to_be_focused()
        expect(page.locator("#config-yaml-lines .yaml-line-refused")).to_contain_text(str(expected_line))
        expect(page.locator("#config-yaml-apply")).to_be_disabled()
        expect(page.locator("#config-yaml-download")).to_be_disabled()
        assert page.evaluate("document.scrollingElement.scrollWidth <= innerWidth")
        page.keyboard.press("Escape")
        expect(page.locator("#config-yaml-dialog")).to_be_hidden()
        expect(page.locator("#config-yaml-open")).to_be_focused()
        assert not page.errors
    finally: page.context.close()


def test_validated_yaml_applies_exact_choices_only_after_the_explicit_action(site):
    page = editor(site)
    try:
        before = page.evaluate("window.FastMDXRun.currentState()")
        text = document(site, 317)
        page.locator("#config-yaml-text").fill(text)
        page.locator("#config-yaml-validate").click()
        expect(page.locator("#config-yaml-apply")).to_be_enabled()
        assert page.evaluate("window.FastMDXRun.currentState()") == before
        with page.expect_download() as downloaded:
            page.locator("#config-yaml-download").click()
        assert open(downloaded.value.path()).read() == text
        page.locator("#config-yaml-apply").click()
        expect(page.locator("#config-yaml-dialog")).to_be_hidden()
        built = page.evaluate("window.FastMDXRun.fetchConfig()")
        assert yaml.safe_load(built["yaml"])["simulation"]["temperature_K"] == 317
        assert not page.errors
    finally: page.context.close()


def test_unrepresentable_valid_settings_are_downloadable_and_never_silently_lost(site):
    page = editor(site)
    try:
        text = document(site) + "analysis:\n  stride: 5\n"
        page.locator("#config-yaml-text").fill(text)
        page.locator("#config-yaml-validate").click()
        expect(page.locator("#config-yaml-status")).to_contain_text("Builder cannot retain")
        expect(page.locator("#config-yaml-apply")).to_be_disabled()
        expect(page.locator("#config-yaml-download")).to_be_enabled()
        assert not page.errors
    finally: page.context.close()


def test_draft_survives_reload_and_form_changes_until_explicit_replacement(site):
    page = editor(site)
    try:
        text = document(site, 333) + "# unsaved form-independent note\n"
        page.locator("#config-yaml-text").fill(text)
        page.reload(wait_until="domcontentloaded")
        page.wait_for_function("window.FastMDXRun.ready")
        page.locator("#config-yaml-open").click()
        expect(page.locator("#config-yaml-text")).to_have_value(text)
        expect(page.locator("#config-yaml-apply")).to_be_disabled()
        page.locator("#config-yaml-replace").click()
        page.locator("#config-yaml-replace-no").click()
        expect(page.locator("#config-yaml-text")).to_have_value(text)
        page.locator("#config-yaml-replace").click()
        page.locator("#config-yaml-replace-yes").click()
        expect(page.locator("#config-yaml-text")).not_to_have_value(text)
        assert not page.errors
    finally: page.context.close()


def test_late_validation_cannot_enable_apply_for_newer_text_or_a_closed_dialog(site):
    page = editor(site)
    try:
        text = document(site, 311)
        page.locator("#config-yaml-text").fill(text)
        answer = page.request.post(site["url"].split("#")[0].rstrip("/") + "/api/config-yaml", data={"text":text}).json()
        held = []
        page.route("**/api/config-yaml", lambda route: held.append(route))
        with page.expect_request("**/api/config-yaml"):
            page.locator("#config-yaml-validate").click()
        page.locator("#config-yaml-text").fill(document(site, 322))
        held.pop().fulfill(content_type="application/json", body=json.dumps(answer))
        expect(page.locator("#config-yaml-apply")).to_be_disabled()
        expect(page.locator("#config-yaml-status")).to_contain_text("changed text")
        with page.expect_request("**/api/config-yaml"):
            page.locator("#config-yaml-validate").click()
        page.keyboard.press("Escape")
        held.pop().fulfill(content_type="application/json", body=json.dumps(answer))
        expect(page.locator("#config-yaml-dialog")).to_be_hidden()
        assert not page.errors
    finally: page.context.close()


def test_browser_mapping_keeps_per_system_choices_and_decision_provenance(site):
    page = editor(site)
    try:
        data = {"systems":[{"system":str(site["root"] / "peptide.pdb"), "id":"one",
                             "simulation":{"temperature_K":318}}],
                "include_phase":["setup","simulation"],"simulation":{"temperature_K":310,"random_seed":42},
                "decisions":{"simulation.temperature_K":{"why":"Keep the requested temperature."}}}
        page.locator("#config-yaml-text").fill(yaml.safe_dump(data))
        page.locator("#config-yaml-validate").click()
        expect(page.locator("#config-yaml-apply")).to_be_enabled()
        page.locator("#config-yaml-apply").click()
        built = page.evaluate("window.FastMDXRun.fetchConfig()")
        actual = yaml.safe_load(built["yaml"])
        assert actual["systems"][0]["simulation"]["temperature_K"] == 318
        assert actual["simulation"]["temperature_K"] == 310
        assert actual["simulation"]["random_seed"] == 42
        assert actual["decisions"] == data["decisions"]
    finally: page.context.close()


def test_storage_failure_holds_the_workspace_switch_and_keeps_typed_text(site):
    page = editor(site)
    try:
        page.evaluate("""() => { const write = Storage.prototype.setItem;
          Storage.prototype.setItem = function(k,v) { if(k.startsWith('fmx.yamlDraft')) throw new Error('test storage full'); return write.call(this,k,v); }; }""")
        text = document(site, 323)
        page.locator("#config-yaml-text").fill(text)
        expect(page.locator("#config-yaml-status")).to_contain_text("could not save")
        answer = page.evaluate("window.FastMDXDashboard.openStudyTab('aaaaaaaaaaaaaaaa','overview')")
        assert not answer["ok"] and "YAML" in answer["error"]
        expect(page.locator("#config-yaml-text")).to_have_value(text)
        page.keyboard.press("Escape")
        page.locator("#config-yaml-open").click()
        expect(page.locator("#config-yaml-text")).to_have_value(text)
    finally: page.context.close()


def test_newer_shared_draft_is_never_overwritten_by_the_old_window(site):
    page = editor(site)
    try:
        page.locator("#config-yaml-text").fill(document(site, 321))
        remote = json.dumps({"text":document(site, 322)})
        page.evaluate("""remote => {const key=Object.keys(localStorage).find(k=>k.startsWith('fmx.yamlDraft'));localStorage.setItem(key,remote);}""", remote)
        text = document(site, 333)
        page.locator("#config-yaml-text").fill(text)
        expect(page.locator("#config-yaml-status")).to_contain_text("Another window changed")
        assert "YAML" in page.evaluate("window.FastMDXRun.beforeStudySwitch()")
        assert page.evaluate("localStorage.getItem(Object.keys(localStorage).find(k=>k.startsWith('fmx.yamlDraft')))") == remote
        expect(page.locator("#config-yaml-text")).to_have_value(text)
    finally: page.context.close()


def test_pending_builder_copy_cannot_replace_text_typed_while_it_waits(site):
    page = editor(site)
    try:
        built = page.evaluate("window.FastMDXRun.fetchConfig()")
        held = []
        page.route("**/api/config", lambda route: held.append(route))
        page.locator("#config-yaml-replace").click()
        with page.expect_request("**/api/config"):
            page.locator("#config-yaml-replace-yes").click()
        text = document(site, 334)
        page.locator("#config-yaml-text").fill(text)
        held.pop().fulfill(content_type="application/json", body=json.dumps(built))
        expect(page.locator("#config-yaml-text")).to_have_value(text)
        expect(page.locator("#config-yaml-apply")).to_be_disabled()
    finally: page.context.close()


@pytest.mark.parametrize("shortcut", ["Control+Enter", "Meta+Enter"])
def test_yaml_dialog_never_launches_the_hidden_builder_with_its_run_shortcut(site, shortcut):
    page = editor(site)
    try:
        expect(page.locator("#run-start-button")).to_be_enabled()
        page.locator("#config-yaml-text").fill("simulation: {temperature_K: 311}\n")
        page.locator("#config-yaml-validate").click()
        expect(page.locator("#config-yaml-status")).to_contain_text("systems")
        starts = []
        def intercepted(route):
            starts.append(route.request.post_data_json)
            route.fulfill(content_type="application/json", body=json.dumps({"ok":False,"error":"Test launch intercepted."}))
        page.route("**/api/run", intercepted)
        page.locator("#config-yaml-text").focus()
        page.keyboard.press(shortcut)
        page.wait_for_timeout(250)
        assert starts == []
        expect(page.locator("#config-yaml-dialog")).to_be_visible()
        expect(page.locator("#config-yaml-apply")).to_be_disabled()
    finally: page.context.close()


@pytest.mark.parametrize("name,value", [("include",["setup","simulation"]),
    ("include_phase","setup, simulation"), ("exclude",["analysis","report"])])
def test_phase_aliases_apply_the_same_selected_phases_in_the_actual_builder(site, name, value):
    page = editor(site)
    try:
        config = yaml.safe_load(document(site, 319))
        config.pop("include_phase")
        config[name] = value
        page.locator("#config-yaml-text").fill(yaml.safe_dump(config))
        page.locator("#config-yaml-validate").click()
        expect(page.locator("#config-yaml-apply")).to_be_enabled()
        page.locator("#config-yaml-apply").click()
        actual = yaml.safe_load(page.evaluate("window.FastMDXRun.fetchConfig()")["yaml"])
        assert actual["include_phase"] == ["setup","simulation"]
        assert actual["simulation"]["temperature_K"] == 319
    finally: page.context.close()
