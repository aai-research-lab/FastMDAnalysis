"""Open study documents keep their places without closing a simulation."""
from __future__ import annotations

import json
import re
import shutil
import subprocess
import sys

import pytest

sync_api = pytest.importorskip("playwright.sync_api")
expect = sync_api.expect


def _study(path, name):
    (path / "simulation").mkdir(parents=True)
    (path / "manifest.json").write_text(json.dumps({"system": name, "phases": [
        {"name": "simulation", "status": "completed"}]}))
    (path / "simulation" / "live_status.json").write_text(json.dumps({
        "status": "completed", "stage": "production", "system_id": name,
        "current_step": 10, "total_steps": 10}))
    return path


@pytest.fixture(scope="module")
def browser():
    with sync_api.sync_playwright() as pw:
        launched = pw.chromium.launch(args=["--enable-unsafe-swiftshader"])
        yield launched
        launched.close()


@pytest.fixture
def site(tmp_path, monkeypatch):
    from fastmdxplora.gui.server import start_dashboard_session
    from fastmdxplora.gui.workspace import study_key

    monkeypatch.setenv("FASTMDXPLORA_CONFIG_DIR", str(tmp_path / "config"))
    root = tmp_path / "workspace"
    studies = [_study(root / "first", "Alpha"), _study(root / "second", "Beta"),
               _study(root / "third", "Gamma")]
    session = start_dashboard_session(output=str(studies[0]), exploration_root=root,
                                      host="127.0.0.1", port=0)
    yield {"session": session, "studies": studies, "keys": [study_key(p) for p in studies], "root": root}
    session.stop()


def _open(browser, site, width=1440, theme="dark", context=None):
    context = context or browser.new_context(viewport={"width": width, "height": 900})
    context.add_init_script(f"localStorage.setItem('fmx.theme', '{theme}');")
    page = context.new_page()
    page.set_default_timeout(20000)
    errors = []
    page.on("pageerror", lambda error: errors.append(str(error)))
    page.errors = errors
    page.goto(site["session"].url + "#overview", wait_until="domcontentloaded")
    page.locator('.study-tab-link[aria-current="page"]').wait_for()
    return page


def _link(page, key):
    return page.locator(f'.study-tab-link[data-study-key="{key}"]')


def _seen(page):
    return page.locator(".study-tab-link").evaluate_all("es => es.map(e => e.dataset.studyKey)")


def _active(page, key):
    # A switch loads the page again: on a busy machine that takes more than
    # the 5 s an expectation waits by default.
    expect(_link(page, key)).to_have_attribute("aria-current", "page", timeout=20000)
    page.wait_for_function("key => new URL(location.href).searchParams.get('study') === key", arg=key,
                           timeout=20000)


def _outside_open(page, site, index, route="overview"):
    page.evaluate("""async ({path, route}) => {
        const answer = await (await fetch('/api/explore/switch', {
            method: 'POST', headers: {'content-type': 'application/json'},
            body: JSON.stringify({folder: path})})).json();
        window.FastMDXDashboard.studyOpened(answer, route);
    }""", {"path": str(site["studies"][index]), "route": route})
    _active(page, site["keys"][index])


def test_tabs_restore_each_studys_last_page_and_survive_reload(browser, site):
    page = _open(browser, site)
    try:
        page.evaluate("window.FastMDXDashboard.navigate('analysis/rmsd')")
        page.wait_for_function("location.hash === '#analysis/rmsd'")
        _outside_open(page, site, 1)
        assert _seen(page) == site["keys"][:2]
        _link(page, site["keys"][0]).click()
        _active(page, site["keys"][0])
        expect(_link(page, site["keys"][0])).to_be_focused()
        assert page.url.endswith("#analysis/rmsd")
        page.reload(wait_until="domcontentloaded")
        _active(page, site["keys"][0])
        assert _seen(page) == site["keys"][:2]
        _link(page, site["keys"][1]).click()
        _active(page, site["keys"][1])
        assert page.url.endswith("#overview")
        assert page.errors == []
    finally:
        page.context.close()


def test_tab_name_updates_when_the_results_supply_the_study_name(browser, site):
    (site["studies"][0] / "simulation" / "live_status.json").unlink()
    context = browser.new_context()
    def results(route):
        payload = route.fetch().json()
        payload["system"] = {"system": "2XYZ"}
        route.fulfill(json=payload)
    context.route("**/api/results", results)
    page = _open(browser, site, context=context)
    try:
        expect(page.locator("#topbar-run-title")).to_have_text("2XYZ")
        expect(_link(page, site["keys"][0])).to_have_text("2XYZ")
        assert page.errors == []
    finally:
        page.context.close()


def test_a_switch_saves_the_page_before_its_navigation_frame_runs(browser, site):
    page = _open(browser, site)
    try:
        _outside_open(page, site, 1)
        page.evaluate("""key => {
            window.FastMDXDashboard.navigate('analysis/rmsd');
            document.querySelector('.study-tab-link[data-study-key="' + key + '"]').click();
        }""", site["keys"][0])
        _active(page, site["keys"][0])
        _link(page, site["keys"][1]).click()
        _active(page, site["keys"][1])
        assert page.url.endswith("#analysis/rmsd")
        assert page.errors == []
    finally:
        page.context.close()


def test_close_inactive_changes_no_study_and_close_active_selects_neighbour(browser, site):
    page = _open(browser, site)
    try:
        _outside_open(page, site, 1)
        _outside_open(page, site, 2)
        requests = []
        page.on("request", lambda request: requests.append(request.url) if '/api/explore/switch' in request.url else None)
        page.locator(f'.study-tab-close[data-study-key="{site["keys"][0]}"]').click()
        assert _seen(page) == site["keys"][1:]
        assert requests == []
        page.locator(f'.study-tab-close[data-study-key="{site["keys"][2]}"]').click()
        _active(page, site["keys"][1])
        assert _seen(page) == [site["keys"][1]]
        assert all((p / "manifest.json").is_file() for p in site["studies"])
        assert page.errors == []
    finally:
        page.context.close()


def test_last_closed_tab_stays_closed_until_explicitly_opened_again(browser, site):
    page = _open(browser, site)
    try:
        page.locator(".study-tab-close").click()
        page.wait_for_function("location.hash === '#studies'")
        assert _seen(page) == []
        page.evaluate("""() => window.FastMDXDashboard.applyStatus({
            status: {system_id:'Alpha',status:'completed'}, study: new URL(location).searchParams.get('study')})""")
        assert _seen(page) == []
        page.reload(wait_until="domcontentloaded")
        page.locator("#study-tabs:not([hidden])").wait_for()
        assert _seen(page) == []
        _outside_open(page, site, 0)
        assert _seen(page) == [site["keys"][0]]
        assert page.errors == []
    finally:
        page.context.close()


def test_a_refused_switch_or_close_keeps_the_current_document(browser, site):
    page = _open(browser, site)
    try:
        _outside_open(page, site, 1)
        before = page.url
        page.route("**/api/explore/switch", lambda route: route.fulfill(
            status=403, content_type="application/json", body=json.dumps({"ok": False, "error": "Switch refused for this check."})))
        _link(page, site["keys"][0]).click()
        expect(page.locator("#study-tabs-note")).to_have_text("Switch refused for this check.")
        assert page.url == before
        _active(page, site["keys"][1])
        page.locator(f'.study-tab-close[data-study-key="{site["keys"][1]}"]').click()
        expect(page.locator("#study-tabs")).to_have_attribute("aria-busy", "false")
        assert _seen(page) == site["keys"][:2]
        assert page.url == before
        assert page.errors == []
    finally:
        page.context.close()


def test_pending_switch_starts_one_request_and_does_not_relabel_old_records(browser, site):
    page = _open(browser, site)
    try:
        _outside_open(page, site, 1)
        held = []
        page.route("**/api/explore/switch", lambda route: held.append(route))
        _link(page, site["keys"][0]).click()
        expect(page.locator("#study-tabs")).to_have_attribute("aria-busy", "true")
        _link(page, site["keys"][1]).dispatch_event("click", {"button": 0})
        assert len(held) == 1
        _active(page, site["keys"][1])
        held[0].fulfill(status=500, content_type="application/json", body=json.dumps({"ok": False, "error": "Try later."}))
        expect(page.locator("#study-tabs-note")).to_have_text("Try later.")
        assert page.errors == []
    finally:
        page.context.close()


def test_browser_history_returns_to_the_study_and_its_page(browser, site):
    page = _open(browser, site)
    try:
        page.evaluate("window.FastMDXDashboard.navigate('analysis/rmsd')")
        _outside_open(page, site, 1, "files")
        _link(page, site["keys"][0]).click()
        _active(page, site["keys"][0])
        assert page.url.endswith("#analysis/rmsd")
        page.go_back()
        _active(page, site["keys"][1])
        assert page.url.endswith("#files")
        page.go_forward()
        _active(page, site["keys"][0])
        assert page.url.endswith("#analysis/rmsd")
        assert _seen(page) == site["keys"][:2]
        assert page.errors == []
    finally:
        page.context.close()


def test_untrusted_stored_entries_are_text_and_only_admitted_keys_are_requested(browser, site):
    context = browser.new_context(viewport={"width": 1440, "height": 900})
    key = "fmx.studyTabs:" + str(site["root"])
    evil = '<img src=x onerror="window.badTab=true">'
    entries = {"tabs": [{"key": site["keys"][1], "title": evil, "page": "https://attacker.invalid/"},
                         {"key": site["keys"][1], "title": "duplicate", "page": "overview"},
                         {"key": "../../private", "title": "Invalid", "page": "overview"}], "closed": []}
    context.add_init_script("sessionStorage.setItem(" + json.dumps(key) + ", " + json.dumps(json.dumps(entries)) + ");")
    page = _open(browser, site, context=context)
    try:
        assert _seen(page) == [site["keys"][1], site["keys"][0]]
        assert _link(page, site["keys"][1]).text_content() == evil
        assert page.locator("#study-tabs img").count() == 0
        assert page.evaluate("!window.badTab")
        assert _link(page, site["keys"][1]).get_attribute("href").endswith("#overview")
        assert page.errors == []
    finally:
        page.context.close()


def test_an_unsent_agent_message_holds_the_switch_without_discarding_it(browser, site):
    page = _open(browser, site)
    try:
        _outside_open(page, site, 1, "agent")
        page.fill("#agent-request", "Keep this unsent question.")
        requests = []
        page.on("request", lambda request: requests.append(request.url) if '/api/explore/switch' in request.url else None)
        _link(page, site["keys"][0]).click()
        expect(page.locator("#study-tabs-note")).to_contain_text("unsent Agent message")
        assert requests == []
        assert page.input_value("#agent-request") == "Keep this unsent question."
        _active(page, site["keys"][1])
        page.fill("#agent-request", "")
        _link(page, site["keys"][0]).click()
        _active(page, site["keys"][0])
        assert page.errors == []
    finally:
        page.context.close()


def test_an_unsent_attachment_holds_the_switch_until_removed(browser, site):
    page = _open(browser, site)
    try:
        _outside_open(page, site, 1, "agent")
        path = str(site["studies"][1] / "manifest.json")
        page.locator("#agent-attach-path").evaluate("""(input, path) => {
            input.value = path; input.dispatchEvent(new Event('change'));
        }""", path)
        page.locator("#agent-attachments .agent-chip").wait_for()
        _link(page, site["keys"][0]).click()
        expect(page.locator("#study-tabs-note")).to_contain_text("attachments")
        _active(page, site["keys"][1])
        expect(page.locator("#agent-attachments .agent-chip")).to_be_visible()
        page.get_by_role("button", name="Remove manifest.json", exact=True).click()
        _link(page, site["keys"][0]).click()
        _active(page, site["keys"][0])
        assert page.errors == []
    finally:
        page.context.close()


@pytest.mark.parametrize("streaming", [True, False])
def test_an_agent_reply_holds_the_switch_until_it_finishes(browser, site, streaming):
    page = _open(browser, site)
    try:
        _outside_open(page, site, 1, "agent")
        if not streaming:
            page.evaluate("window.ReadableStream = undefined")
        held = []
        page.route("**/api/agent/propose*", lambda route: held.append(route))
        page.fill("#agent-request", "Describe this study.")
        page.keyboard.press("Enter")
        page.locator("#agent-propose.is-writing").wait_for()
        _link(page, site["keys"][0]).click()
        expect(page.locator("#study-tabs-note")).to_contain_text("Agent reply")
        _active(page, site["keys"][1])
        assert len(held) == 1
        held[0].fulfill(json={"ok": False, "answer": "The reply has finished.", "cites": []})
        page.get_by_text("The reply has finished.", exact=True).wait_for()
        _link(page, site["keys"][0]).click()
        _active(page, site["keys"][0])
        assert page.errors == []
    finally:
        page.context.close()


def test_a_folder_without_study_records_never_becomes_a_study_tab(browser, site):
    from fastmdxplora.gui.server import start_dashboard_session

    empty = site["root"] / "empty"
    empty.mkdir()
    session = start_dashboard_session(output=str(empty), exploration_root=site["root"],
                                      host="127.0.0.1", port=0)
    context = browser.new_context()
    page = context.new_page()
    try:
        page.goto(session.url + "#overview", wait_until="domcontentloaded")
        page.wait_for_function("""() => {
            const state = window.FastMDXDashboard?.state;
            return state?.studyIsOne === false && state.results?.output_dir;
        }""")
        expect(page.locator("#study-tabs")).to_be_visible()
        assert _seen(page) == []
        page.evaluate("window.FastMDXDashboard.navigate('files')")
        page.wait_for_function("location.hash === '#files'")
        # Navigation is announced on the next frame.
        page.evaluate("() => new Promise(requestAnimationFrame)")
        assert _seen(page) == []
        assert page.locator('.study-tab-link[aria-current="page"]').count() == 0
    finally:
        context.close()
        session.stop()


def test_an_attachment_still_loading_holds_the_switch(browser, site):
    page = _open(browser, site)
    try:
        _outside_open(page, site, 1, "agent")
        held = []
        page.route("**/api/agent/attachment", lambda route: held.append(route))
        with page.expect_request("**/api/agent/attachment"):
            page.locator("#agent-attach-path").evaluate("""input => {
                input.value = '/pending/notes.txt'; input.dispatchEvent(new Event('change'));
            }""")
        _link(page, site["keys"][0]).click()
        expect(page.locator("#study-tabs-note")).to_contain_text("attachment to finish loading")
        _active(page, site["keys"][1])
        held[0].fulfill(json={"ok": True, "path": "/pending/notes.txt", "name": "notes.txt",
                              "size": 4, "text": "kept"})
        expect(page.locator("#agent-attachments .agent-chip")).to_contain_text("notes.txt")
        assert page.errors == []
    finally:
        page.context.close()


def _builder(page, site):
    pdb = site["root"] / "peptide.pdb"
    pdb.write_text("ATOM      1  N   GLY A   1       0.000   0.000   0.000  1.00  0.00           N\n"
                   "ATOM      2  CA  GLY A   1       1.450   0.000   0.000  1.00  0.00           C\n"
                   "ATOM      3  C   GLY A   1       2.900   0.000   0.000  1.00  0.00           C\nEND\n")
    page.evaluate("window.FastMDXDashboard.navigate('run')")
    page.click("#run-start [data-start='structure']")
    page.fill("#run-system", str(pdb))
    page.locator('#run-settings [data-setting="ph"] input').fill("6.5")


def test_a_builder_launch_waiting_for_its_answer_holds_the_switch(browser, site):
    page = _open(browser, site)
    try:
        _outside_open(page, site, 1)
        _builder(page, site)
        held = []
        page.route("**/api/run", lambda route: held.append(route))
        with page.expect_request("**/api/run"):
            page.click("#run-start-button")
        _link(page, site["keys"][0]).click()
        expect(page.locator("#study-tabs-note")).to_contain_text("Config launch to finish")
        _active(page, site["keys"][1])
        held[0].fulfill(json={"ok": False, "error": "Launch refused for this check."})
        expect(page.locator("#run-note")).to_have_text("Launch refused for this check.")
        _link(page, site["keys"][0]).click()
        _active(page, site["keys"][0])
        assert page.errors == []
    finally:
        page.context.close()


def test_switch_flushes_and_restores_the_workspace_builder_draft(browser, site):
    page = _open(browser, site)
    try:
        _outside_open(page, site, 1)
        _builder(page, site)
        _link(page, site["keys"][0]).click()
        _active(page, site["keys"][0])
        page.evaluate("window.FastMDXDashboard.navigate('run')")
        expect(page.locator('#run-settings [data-setting="ph"] input')).to_have_value("6.5")
        assert page.errors == []
    finally:
        page.context.close()


def test_a_failed_draft_save_keeps_the_current_study_and_typed_value(browser, site):
    page = _open(browser, site)
    try:
        _outside_open(page, site, 1)
        _builder(page, site)
        page.evaluate("""() => {
            const original = Storage.prototype.setItem;
            Storage.prototype.setItem = function(k,v) {
                if (k.startsWith('fmx.builderDraft')) throw new DOMException('full','QuotaExceededError');
                return original.call(this,k,v);
            };
        }""")
        page.locator('#run-settings [data-setting="ph"] input').fill("6.6")
        _link(page, site["keys"][0]).click()
        expect(page.locator("#study-tabs-note")).to_contain_text("draft could not be saved")
        _active(page, site["keys"][1])
        assert page.input_value('#run-settings [data-setting="ph"] input') == "6.6"
        assert page.errors == []
    finally:
        page.context.close()


LAYOUT = """() => {
    const box = el => el && el.getClientRects().length ? el.getBoundingClientRect() : null;
    const visible = el => el && el.getClientRects().length && getComputedStyle(el).visibility !== 'hidden';
    const host = document.getElementById('study-tabs'), list = document.getElementById('study-tabs-list');
    const strip = box(host), room = box(list), current = box(list.querySelector('.study-tab.is-current'));
    const open = box(document.getElementById('study-tabs-open'));
    const panel = [...document.querySelectorAll('.side-panel, .agent-drawer')].filter(visible).map(box)
        .filter(r => r.width && r.left > innerWidth / 2);
    const title = box(document.querySelector('.page.active .page-title, .page[data-page="overview"] .page-title'));
    return {
        wider: document.documentElement.scrollWidth - innerWidth,
        strip: [strip.left, strip.right, strip.top, strip.bottom],
        currentInView: Boolean(current) && current.left >= room.left - 1 && current.right <= room.right + 1,
        openAfterList: room.right <= open.left + 1 && open.right <= strip.right + 1,
        scrolls: list.scrollWidth > list.clientWidth,
        moreShown: host.classList.contains('more-before') || host.classList.contains('more-after'),
        panelLeft: panel.length ? Math.min(...panel.map(r => r.left)) : innerWidth,
        titleTop: title ? title.top : null,
        stops: [...host.querySelectorAll('a, button')].filter(e => e.tabIndex >= 0).map(e => e.className),
    };
}"""


@pytest.mark.parametrize("width,beside", [(1440, False), (1250, True), (900, False), (390, False)])
@pytest.mark.parametrize("theme", ["light", "dark"])
def test_many_tabs_fit_the_page_and_work_by_keyboard(browser, site, width, beside, theme):
    context = browser.new_context(viewport={"width": width, "height": 900})
    key = "fmx.studyTabs:" + str(site["root"])
    entries = {"tabs": [{"key": site["keys"][1], "title": "Beta" * 60, "page": "overview"}]
               + [{"key": f"{n:016x}", "title": f"Study {n}", "page": "overview"} for n in range(10, 28)]
               + [{"key": site["keys"][0], "title": "Alpha", "page": "overview"}]}
    context.add_init_script("sessionStorage.setItem(" + json.dumps(key) + ", " + json.dumps(json.dumps(entries)) + ");")
    page = _open(browser, site, width=width, theme=theme, context=context)
    try:
        if beside:
            page.click("#agent-beside-open")
            page.locator("#agent-drawer").wait_for()
        page.wait_for_function("document.getElementById('study-tabs').matches('.more-before, .more-after')")
        facts = page.evaluate(LAYOUT)
        # The open study is in view at the end of a long strip; a fade says
        # more tabs lie past the edge; nothing runs past the page or under
        # the side panel or the Agent, and the page's title starts below.
        assert facts["wider"] <= 0, facts
        assert facts["strip"][0] >= 0 and facts["strip"][1] <= min(width, facts["panelLeft"]) + 1, facts
        assert facts["currentInView"] and facts["scrolls"] and facts["moreShown"], facts
        assert facts["openAfterList"], facts
        assert facts["titleTop"] is None or facts["titleTop"] >= facts["strip"][3] - 1, facts
        # One keyboard stop for the strip's tabs, then its button.
        assert facts["stops"] == ["study-tab-link", "line-btn study-tabs-open"], facts
        expect(_link(page, site["keys"][0])).to_have_attribute("tabindex", "0")
        _link(page, site["keys"][0]).focus()
        page.keyboard.press("Home")
        expect(_link(page, site["keys"][1])).to_be_focused()
        page.keyboard.press("End")
        expect(_link(page, site["keys"][0])).to_be_focused()
        page.keyboard.press("Home")
        # The current tab and the focused one look different.
        looks = page.evaluate("""keys => keys.map(k => {
            const item = document.querySelector('.study-tab-link[data-study-key="' + k + '"]').parentElement;
            const css = getComputedStyle(item);
            return css.backgroundColor + css.boxShadow;
        })""", [site["keys"][0], site["keys"][1]])
        assert looks[0] != looks[1]
        page.keyboard.press("Delete")
        assert site["keys"][1] not in _seen(page)
        assert page.evaluate("document.documentElement.scrollWidth <= innerWidth")
        assert page.errors == []
    finally:
        page.context.close()


def test_a_delayed_success_does_not_show_new_records_under_the_old_tab(browser, site):
    page = _open(browser, site)
    try:
        _outside_open(page, site, 1)
        held = []
        def switch(route):
            response = route.fetch()
            held.append((route, response))
        page.route("**/api/explore/switch", switch)
        _link(page, site["keys"][0]).click()
        expect(page.locator("#study-tabs")).to_have_attribute("aria-busy", "true")
        assert held
        assert page.evaluate("getComputedStyle(document.querySelector('.page-shell')).visibility") == "hidden"
        assert page.evaluate("window.FastMDXDashboard.state.studyShown") == site["keys"][1]
        # The actual backend has changed while the HTTP reply is held.
        server_key = page.evaluate("async () => (await (await fetch('/api/status')).json()).study")
        assert server_key == site["keys"][0]
        _active(page, site["keys"][1])
        held[0][0].fulfill(response=held[0][1])
        _active(page, site["keys"][0])
        expect(page.locator(".page-shell")).to_have_css("visibility", "visible")
        assert page.errors == []
    finally:
        page.context.close()


def _agent_versions(page):
    import yaml

    from fastmdxplora.gui.plan import changes_between

    base = {"systems": [{"system": "1UBQ"}], "include_phase": ["simulation"],
            "simulation": {"temperature_K": 300}}
    after = {**base, "simulation": {"temperature_K": 330}}
    entries = [{"role": "agent", "kind": "config", "config": c,
                "yaml": yaml.safe_dump(c), "plan": [], "attempts": [], "cycles": 1,
                "version": i + 1} for i, c in enumerate([base, after])]
    entries[1]["changes"] = changes_between(base, after)
    page.evaluate("""async entries => {
        await fetch('/api/agent/conversation', {method:'POST',
            headers:{'content-type':'application/json'}, body:JSON.stringify({entries})});
    }""", entries)
    page.reload(wait_until="domcontentloaded")
    return page.locator('.agent-study:not(.is-replaced):not([hidden])')


def test_a_rejected_change_not_applied_holds_the_tab_switch(browser, site):
    page = _open(browser, site)
    try:
        _outside_open(page, site, 1, "agent")
        current = _agent_versions(page)
        row = current.locator('[data-review-setting="simulation.temperature_K"]')
        row.locator('[data-choice="reject"]').click()
        requests = []
        page.on("request", lambda request: requests.append(request.url) if '/api/explore/switch' in request.url else None)
        _link(page, site["keys"][0]).click()
        # The reason names the review's own buttons, here and beside it.
        said = "Use reviewed version, or Use this version to keep version 2 as proposed, before switching studies."
        expect(page.locator("#study-tabs-note")).to_have_text(said)
        expect(current.locator(".agent-review-note")).to_contain_text("before switching studies")
        assert requests == []
        _active(page, site["keys"][1])
        assert row.locator('[data-choice="reject"]').get_attribute("aria-pressed") == "true"
        current.get_by_role("button", name="Use this version", exact=True).click()
        page.locator('.agent-study:not(.is-replaced):not([hidden]) [data-role="version"]').filter(has_text="3").wait_for()
        _link(page, site["keys"][0]).click()
        _active(page, site["keys"][0])
        assert page.errors == []
    finally:
        page.context.close()


def test_an_agent_launch_waiting_for_its_answer_holds_the_switch(browser, site):
    page = _open(browser, site)
    try:
        _outside_open(page, site, 1, "agent")
        current = _agent_versions(page)
        held = []
        page.route("**/api/agent/run", lambda route: held.append(route))
        with page.expect_request("**/api/agent/run"):
            current.locator('[data-role="run"]').click()
        _link(page, site["keys"][0]).click()
        expect(page.locator("#study-tabs-note")).to_contain_text("study launch to finish")
        _active(page, site["keys"][1])
        held[0].fulfill(json={"ok": False, "error": "Launch refused for this check."})
        page.get_by_text("Launch refused for this check.", exact=True).first.wait_for()
        _link(page, site["keys"][0]).click()
        _active(page, site["keys"][0])
        assert page.errors == []
    finally:
        page.context.close()


def _switches(page):
    asked = []
    page.on("request", lambda request: asked.append(request.url)
            if request.method == "POST" and "/api/explore/switch" in request.url else None)
    return asked


def test_a_gone_folder_of_another_tab_is_said_on_it_and_not_opened(browser, site):
    page = _open(browser, site)
    try:
        _outside_open(page, site, 1)
        before = page.url
        shutil.rmtree(site["studies"][0])
        _link(page, site["keys"][0]).click()
        note = page.locator("#study-tabs-note")
        expect(note).to_have_text("This study's folder is no longer there.")
        # No path of the server's is said.
        assert str(site["root"]) not in note.text_content()
        expect(_link(page, site["keys"][0])).to_contain_text("folder gone")
        assert page.url == before
        _active(page, site["keys"][1])
        # Said on its tab after a reload too; choosing it again asks nothing.
        page.reload(wait_until="domcontentloaded")
        _active(page, site["keys"][1])
        expect(_link(page, site["keys"][0])).to_contain_text("folder gone")
        asked = _switches(page)
        _link(page, site["keys"][0]).click()
        expect(note).to_have_text("This study's folder is no longer there.")
        assert asked == [] and page.url == before
        page.locator(f'.study-tab-close[data-study-key="{site["keys"][0]}"]').click()
        assert _seen(page) == [site["keys"][1]]
        assert page.errors == []
    finally:
        page.context.close()


def test_the_open_studys_gone_folder_is_said_and_not_shown_as_an_empty_study(browser, site):
    page = _open(browser, site)
    try:
        _outside_open(page, site, 1)
        shutil.rmtree(site["studies"][1])
        gone = page.locator("#study-gone")
        expect(gone).to_be_visible(timeout=20000)
        expect(gone).to_contain_text("This study's folder is no longer there")
        expect(page.locator(".page-shell")).to_be_hidden()
        expect(_link(page, site["keys"][1])).to_contain_text("folder gone")
        expect(_link(page, site["keys"][1])).to_contain_text(re.compile("beta", re.I))
        # The address still names it, after a reload as well.
        page.reload(wait_until="domcontentloaded")
        _active(page, site["keys"][1])
        expect(gone).to_be_visible()
        expect(page.locator(".page-shell")).to_be_hidden()
        assert str(site["root"]) not in page.locator("#study-tabs").inner_text()
        # All studies is still a page of its own; Close opens the neighbour.
        page.evaluate("window.FastMDXDashboard.navigate('studies')")
        expect(gone).to_be_hidden()
        page.evaluate("window.FastMDXDashboard.navigate('overview')")
        page.locator("#study-gone-close").click()
        _active(page, site["keys"][0])
        expect(gone).to_be_hidden()
        assert _seen(page) == [site["keys"][0]]
        assert page.errors == []
    finally:
        page.context.close()


def test_closing_the_tab_of_a_running_study_leaves_its_run_going(browser, site):
    first = site["studies"][0]
    (first / "simulation" / "live_status.json").write_text(json.dumps({
        "status": "running", "stage": "production", "system_id": "Alpha",
        "current_step": 10, "total_steps": 5000}))
    runtime = site["session"].runtime
    runtime.process = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(120)"],
                                       start_new_session=True)
    runtime.running_root = first
    page = _open(browser, site)
    try:
        posted = []
        page.on("request", lambda request: posted.append(request.url) if request.method == "POST" else None)
        _outside_open(page, site, 1)
        _outside_open(page, site, 0)
        posted.clear()
        # The open one, a running study: its neighbour opens; the run goes on.
        page.locator(f'.study-tab-close[data-study-key="{site["keys"][0]}"]').click()
        _active(page, site["keys"][1])
        # Only the switch to the neighbour: nothing stops or changes the run.
        assert [u for u in posted if "/api/explore/switch" in u] == posted[:1], posted
        assert not any("stop" in url or "/api/run" in url for url in posted), posted
        assert runtime.process.poll() is None
        state = page.evaluate("async () => (await (await fetch('/api/app-state')).json())")
        assert state["running_elsewhere"] == str(first.resolve())
        # Not the open one: nothing is asked of the server at all.
        _outside_open(page, site, 0)
        _link(page, site["keys"][1]).click()
        _active(page, site["keys"][1])
        posted.clear()
        page.locator(f'.study-tab-close[data-study-key="{site["keys"][0]}"]').click()
        assert _seen(page) == [site["keys"][1]]
        assert not any("switch" in url or "stop" in url for url in posted), posted
        assert runtime.process.poll() is None
        assert json.loads((first / "simulation" / "live_status.json").read_text())["status"] == "running"
        assert page.errors == []
    finally:
        page.context.close()
        runtime.process.kill()
        runtime.process.wait()


def _typed(page, site):
    _outside_open(page, site, 1, "agent")
    page.fill("#agent-request", "Unsent words.")
    return _switches(page)


@pytest.mark.parametrize("way", ["recent", "all studies", "load study", "history"])
def test_every_way_to_another_study_is_held_by_unsent_agent_words(browser, site, way):
    page = _open(browser, site)
    try:
        asked = _typed(page, site)
        before = page.url
        if way == "recent":
            page.locator('.sidebar-recent-item:not([aria-current])').first.click()
        elif way == "all studies":
            page.evaluate("window.FastMDXDashboard.navigate('studies')")
            page.locator('.page[data-page="studies"] .study-open').first.click()
            before = page.url
        elif way == "load study":
            page.locator("#load-study-path").evaluate("""(input, path) => {
                input.value = path; input.dispatchEvent(new Event('change'));
            }""", str(site["studies"][0]))
        else:
            # A step of this page naming another study, gone back to.
            page.evaluate("""key => {
                const here = location.href, there = new URL(here);
                there.searchParams.set('study', key);
                history.pushState(null, '', there.href);
                history.pushState(null, '', here);
                history.back();
            }""", site["keys"][0])
            page.wait_for_function("key => new URL(location.href).searchParams.get('study') === key",
                                   arg=site["keys"][1])
        expect(page.locator("#study-tabs-note")).to_contain_text("unsent Agent message")
        assert asked == []
        assert page.input_value("#agent-request") == "Unsent words."
        _active(page, site["keys"][1])
        if way != "history":
            assert page.url == before
        assert page.errors == []
    finally:
        page.context.close()


def test_back_and_forward_between_studies_keep_unsent_agent_words(browser, site):
    page = _open(browser, site)
    try:
        _typed(page, site)
        page.go_back(wait_until="domcontentloaded")
        _active(page, site["keys"][0])
        assert page.input_value("#agent-request") == ""
        page.go_forward(wait_until="domcontentloaded")
        _active(page, site["keys"][1])
        expect(page.locator("#agent-request")).to_have_value("Unsent words.")
        assert page.errors == []
    finally:
        page.context.close()


def _ph(page):
    return page.locator('#run-settings [data-setting="ph"] input')


def _stored_draft(page):
    return page.evaluate("""() => Object.keys(localStorage).filter(k => k.startsWith('fmx.builderDraft'))
        .map(k => localStorage.getItem(k)).join('')""")


def test_a_newer_draft_from_another_window_is_kept_and_does_not_hold_the_switch(browser, site):
    page = _open(browser, site)
    other = page.context.new_page()
    other.set_default_timeout(20000)
    try:
        _outside_open(page, site, 1)
        _builder(page, site)
        # Typed and left: nothing in this window is being typed now.
        page.evaluate("document.activeElement.blur()")
        page.wait_for_function("() => Object.keys(localStorage).some(k => k.startsWith('fmx.builderDraft') && localStorage.getItem(k).includes('6.5'))")
        other.goto(site["session"].url + "#run", wait_until="domcontentloaded")
        expect(_ph(other)).to_have_value("6.5")
        _ph(other).fill("7.2")
        other.keyboard.press("Tab")
        other.wait_for_function("() => Object.keys(localStorage).some(k => k.startsWith('fmx.builderDraft') && localStorage.getItem(k).includes('7.2'))")
        # This window changed nothing since: it takes the newer draft and
        # switches, and the newer draft is the one kept.
        expect(_ph(page)).to_have_value("7.2")
        _link(page, site["keys"][0]).click()
        _active(page, site["keys"][0])
        assert "7.2" in _stored_draft(page)
        # Changes here not saved, a newer draft from the other window: held,
        # and said as that, not as a draft that could not be saved.
        page.evaluate("window.FastMDXDashboard.navigate('run')")
        expect(_ph(page)).to_have_value("7.2")
        _ph(page).fill("6.1")
        other.evaluate("""() => {
            const key = Object.keys(localStorage).find(k => k.startsWith('fmx.builderDraft'));
            localStorage.setItem(key, localStorage.getItem(key).replace('7.2', '7.4'));
        }""")
        _link(page, site["keys"][1]).click()
        note = page.locator("#study-tabs-note")
        expect(note).to_contain_text("another window saved a newer Config draft")
        assert "could not be saved" not in note.text_content()
        _active(page, site["keys"][0])
        assert page.input_value('#run-settings [data-setting="ph"] input') == "6.1"
        assert "7.4" in _stored_draft(page)
        assert page.errors == []
    finally:
        page.context.close()


def test_without_session_storage_the_strip_says_open_studies_are_not_kept(browser, site):
    context = browser.new_context(viewport={"width": 1440, "height": 900})
    context.add_init_script("""Object.defineProperty(window, 'sessionStorage', {
        get() { throw new DOMException('denied', 'SecurityError'); }});""")
    page = _open(browser, site, context=context)
    try:
        expect(page.locator("#study-tabs-unkept")).to_contain_text("Open studies are not kept in this window")
        assert page.errors == []
    finally:
        page.context.close()


def test_tabs_of_one_system_are_named_by_their_folders_too(browser, site):
    for study in site["studies"][1:]:
        (study / "simulation" / "live_status.json").write_text(json.dumps({
            "status": "completed", "stage": "production", "system_id": "Alpha",
            "current_step": 10, "total_steps": 10}))
    page = _open(browser, site)
    try:
        _outside_open(page, site, 1)
        _outside_open(page, site, 2)
        expect(_link(page, site["keys"][0])).to_have_text("ALPH \u00b7 first")
        expect(_link(page, site["keys"][1])).to_have_text("ALPH \u00b7 second")
        expect(_link(page, site["keys"][2])).to_have_text("ALPH \u00b7 third")
        page.locator(f'.study-tab-close[data-study-key="{site["keys"][1]}"]').click()
        page.locator(f'.study-tab-close[data-study-key="{site["keys"][0]}"]').click()
        # Alone with its system, its system names it.
        expect(_link(page, site["keys"][2])).to_have_text("ALPH")
        assert page.errors == []
    finally:
        page.context.close()


def test_closing_the_last_tab_leaves_no_closed_study_named(browser, site):
    page = _open(browser, site)
    try:
        page.locator(".study-tab-close").click()
        page.wait_for_function("location.hash === '#studies'")
        assert "study=" not in page.url
        page.wait_for_timeout(500)
        page.evaluate("() => window.FastMDXDashboard.applyStatus({status: {system_id: 'Alpha', status: 'completed'}, study: window.FastMDXDashboard.state.studyShown})")
        assert "study=" not in page.url
        page.reload(wait_until="domcontentloaded")
        page.locator("#study-tabs:not([hidden])").wait_for()
        page.wait_for_function("window.FastMDXDashboard.state.addressKnown")
        assert "study=" not in page.url and _seen(page) == []
        assert page.errors == []
    finally:
        page.context.close()


def test_a_hold_note_goes_once_another_page_is_shown(browser, site):
    page = _open(browser, site)
    try:
        _typed(page, site)
        _link(page, site["keys"][0]).click()
        note = page.locator("#study-tabs-note")
        expect(note).to_contain_text("unsent Agent message")
        page.evaluate("window.FastMDXDashboard.navigate('overview')")
        expect(note).to_be_hidden()
        assert page.errors == []
    finally:
        page.context.close()
