"""Open study documents keep their places without closing a simulation."""
from __future__ import annotations

import json

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
    page.set_default_timeout(10000)
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
    expect(_link(page, key)).to_have_attribute("aria-current", "page")
    page.wait_for_function("key => new URL(location.href).searchParams.get('study') === key", arg=key)


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


@pytest.mark.parametrize("width", [390, 980, 1600])
@pytest.mark.parametrize("theme", ["light", "dark"])
def test_many_tabs_and_long_names_fit_and_work_by_keyboard(browser, site, width, theme):
    context = browser.new_context(viewport={"width": width, "height": 900})
    key = "fmx.studyTabs:" + str(site["root"])
    entries = {"tabs": [{"key": site["keys"][1], "title": "Beta" * 60, "page": "overview"}]
               + [{"key": f"{n:016x}", "title": f"Study {n}", "page": "overview"} for n in range(10, 28)]}
    context.add_init_script("sessionStorage.setItem(" + json.dumps(key) + ", " + json.dumps(json.dumps(entries)) + ");")
    page = _open(browser, site, width=width, theme=theme, context=context)
    try:
        _link(page, site["keys"][1]).focus()
        page.keyboard.press("End")
        expect(_link(page, site["keys"][0])).to_be_focused()
        page.keyboard.press("Home")
        expect(_link(page, site["keys"][1])).to_be_focused()
        page.keyboard.press("Delete")
        assert site["keys"][1] not in _seen(page)
        assert page.evaluate("document.documentElement.scrollWidth <= innerWidth")
        bounds = page.locator("#study-tabs").evaluate("""host => {
            const box = host.getBoundingClientRect();
            const list = host.querySelector('ul').getBoundingClientRect();
            const open = host.querySelector('#study-tabs-open').getBoundingClientRect();
            return {left: box.left, right: box.right, end: open.right,
                    separated: list.right <= open.left, scrolls: host.querySelector('ul').scrollWidth > list.width};
        }""")
        assert bounds["right"] <= width + 1 and bounds["left"] >= 0
        assert bounds["end"] <= bounds["right"] and bounds["separated"] and bounds["scrolls"]
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

    base = {"systems": [{"system": "1UBQ"}], "include_phase": ["simulation"],
            "simulation": {"temperature_K": 300}}
    after = {**base, "simulation": {"temperature_K": 330}}
    entries = [{"role": "agent", "kind": "config", "config": c,
                "yaml": yaml.safe_dump(c), "plan": [], "attempts": [], "cycles": 1,
                "version": i + 1} for i, c in enumerate([base, after])]
    page.evaluate("""async entries => {
        await fetch('/api/agent/conversation', {method:'POST',
            headers:{'content-type':'application/json'}, body:JSON.stringify({entries})});
    }""", entries)
    page.reload(wait_until="domcontentloaded")
    return page.locator('.agent-study:not(.is-replaced):not([hidden])')


def test_a_pending_config_review_holds_the_tab_switch(browser, site):
    page = _open(browser, site)
    try:
        _outside_open(page, site, 1, "agent")
        current = _agent_versions(page)
        current.get_by_role("button", name="Review each change", exact=True).click()
        current.locator('[data-review-setting="simulation.temperature_K"]').wait_for()
        _link(page, site["keys"][0]).click()
        expect(page.locator("#study-tabs-note")).to_contain_text("Apply or cancel the Config review")
        _active(page, site["keys"][1])
        assert current.locator('[data-review-setting="simulation.temperature_K"]').is_visible()
        current.get_by_role("button", name="Cancel review", exact=True).click()
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
