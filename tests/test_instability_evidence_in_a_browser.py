"""A frozen full-system observation, real atom mapping and a reviewed Config."""
from __future__ import annotations

import io
import json
import time
from pathlib import Path

import pytest

from tests.test_instability_evidence_is_retained import observed as _observed

observed = _observed
sync_api = pytest.importorskip("playwright.sync_api")
expect = sync_api.expect


def _visible_highlight(page):
    from PIL import Image
    import numpy as np

    deadline = time.monotonic() + 8
    while time.monotonic() < deadline:
        pixels = np.asarray(Image.open(io.BytesIO(page.locator("#instability-canvas canvas").first.screenshot())).convert("RGB")).astype(float)
        height, width = pixels.shape[:2]
        # Focus puts the atom centrally; exclude axes and surrounding controls.
        center = pixels[height // 4:3 * height // 4, width // 4:3 * width // 4]
        red, green, blue = (center[:, :, i] for i in range(3))
        if ((green > 100) & (green > red * 1.35) & (green > blue * 1.15)).sum() >= 12:
            return
        page.wait_for_timeout(100)
    pytest.fail("The selected atom has no visible green highlight in the rendered molecule.")


@pytest.fixture(scope="module")
def browser():
    with sync_api.sync_playwright() as pw:
        instance = pw.chromium.launch(args=["--enable-unsafe-swiftshader"])
        yield instance
        instance.close()


@pytest.fixture
def site(observed, monkeypatch, tmp_path):
    from fastmdxplora.gui.server import start_dashboard_session
    from fastmdxplora.gui.workspace import study_key

    root = observed["root"]
    (root / "simulation/live_status.json").write_text(json.dumps({"status": "failed", "stage": "production", "latest_error": "Particle coordinate is NaN"}))
    monkeypatch.setenv("FASTMDXPLORA_CONFIG_DIR", str(tmp_path / "config"))
    session = start_dashboard_session(output=str(root), exploration_root=tmp_path, host="127.0.0.1", port=0)
    calls = []
    monkeypatch.setattr(session.runtime, "_spawn", lambda *a, **kw: calls.append(a))
    yield {"session": session, "root": root, "key": study_key(root), "calls": calls}
    session.stop()


def _open(browser, site, width=1440, theme="dark"):
    context = browser.new_context(viewport={"width": width, "height": 900}, accept_downloads=True)
    context.add_init_script(f"localStorage.setItem('fmx.theme', '{theme}')")
    page = context.new_page(); page.set_default_timeout(15000)
    page.errors = []; page.on("pageerror", lambda error: page.errors.append(str(error)))
    page.goto(site["session"].url + "#viewer", wait_until="domcontentloaded")
    page.locator("#instability-open").click()
    return page


@pytest.mark.parametrize("width,theme", [(390,"dark"),(390,"light"),(980,"dark"),(980,"light"),(1600,"dark"),(1600,"light")])
def test_the_ranked_water_atom_is_highlighted_at_its_retained_coordinates(browser, site, width, theme):
    page = _open(browser, site, width, theme)
    try:
        page.locator("#instability-canvas[data-ready='true']").wait_for()
        expect(page.locator("#instability-ranking button")).to_have_count(4)
        expect(page.locator("#instability-status")).to_contain_text("do not establish the cause")
        page.locator("#instability-ranking button").first.click()
        page.wait_for_function("window.FastMDXInstability.engine?.scene.selected?.length === 1")
        expect(page.locator("#instability-atom-note")).to_contain_text("Highlighted")
        _visible_highlight(page)
        assert page.evaluate("""() => {const e=window.FastMDXInstability.engine, ref=e.components.selected.ref;
          return [...e.plugin.state.data.cells.values()].some(c => c.transform.parent === ref && c.obj?.data?.repr?.props?.ignoreHydrogens === false);}""")
        atom = page.evaluate("() => {const e=window.FastMDXInstability.engine;return e.atoms(e.scene.selected)[0];}")
        assert (atom["chain"], atom["resi"], atom["resn"], atom["atom"]) == ("W", 1, "HOH", "H1")
        assert [atom["x"], atom["y"], atom["z"]] == pytest.approx([4, 2, 1], abs=.005)
        page.locator("#instability-ranking button").last.click()
        page.wait_for_function("window.FastMDXInstability.engine?.atoms(window.FastMDXInstability.engine.scene.selected)[0]?.atom === 'CA'")
        assert page.evaluate("window.FastMDXInstability.engine.atoms(window.FastMDXInstability.engine.scene.selected)[0].icode") == "B"
        assert page.evaluate("document.scrollingElement.scrollWidth <= innerWidth")
        assert page.locator(".instability-box").bounding_box()["width"] <= width
        page.keyboard.press("Escape")
        expect(page.locator("#instability-dialog")).to_be_hidden()
        expect(page.locator("#instability-open")).to_be_focused()
        assert not site["calls"] and page.errors == []
    finally: page.context.close()


def test_review_and_download_preserve_failed_run_and_never_launch(browser, site):
    page = _open(browser, site)
    before = {p.relative_to(site["root"]): p.read_bytes() for p in site["root"].rglob("*") if p.is_file()}
    try:
        page.locator("#instability-canvas[data-ready='true']").wait_for()
        page.locator("#instability-propose").click()
        expect(page.locator("#instability-diff")).to_contain_text("simulation.timestep_fs: 2 → 1")
        expect(page.locator("#instability-proposal-note")).to_contain_text("does not prove")
        with page.expect_download() as downloaded:
            page.locator("#instability-download").click()
        content = Path(downloaded.value.path()).read_text()
        import yaml
        config = yaml.safe_load(content)
        assert config["simulation"]["timestep_fs"] == 1
        assert config["decisions"]["setup.ph"]["source"] == "person"
        assert all((site["root"] / name).read_bytes() == value for name, value in before.items())
        assert not site["calls"]
        page.reload(wait_until="domcontentloaded")
        assert not site["calls"]
        expect(page.locator("#instability-dialog")).to_be_hidden()
    finally: page.context.close()


def test_changed_config_and_snapshot_refuse_without_launching(browser, site):
    page = _open(browser, site)
    try:
        page.locator("#instability-canvas[data-ready='true']").wait_for()
        path = site["root"] / "resolved_config.yml"; path.write_text(path.read_text() + "\n# changed\n")
        page.locator("#instability-propose").click()
        expect(page.locator("#instability-proposal-note")).to_contain_text("changed since")
        expect(page.locator("#instability-download")).to_be_disabled()
        record = json.loads((site["root"] / "simulation/instability.json").read_text())
        snapshot = site["root"] / "simulation" / record["snapshot"]
        snapshot.write_text(snapshot.read_text() + "REMARK changed\n")
        page.locator("#instability-refresh").click()
        expect(page.locator("#instability-status")).to_contain_text("snapshot changed")
        expect(page.locator("#instability-ranking button")).to_have_count(0)
        assert not site["calls"]
    finally: page.context.close()


def test_legacy_run_explicitly_reports_missing_force_evidence(browser, site):
    (site["root"] / "simulation/instability.json").unlink()
    page = _open(browser, site)
    try:
        expect(page.locator("#instability-status")).to_contain_text("Older runs did not record")
        expect(page.locator("#instability-propose")).to_be_disabled()
        expect(page.locator("#instability-canvas")).to_be_hidden()
    finally: page.context.close()


@pytest.mark.parametrize("width", [390, 1440])
def test_problems_opens_the_same_frozen_evidence(browser, site, width):
    page = _open(browser, site, width)
    try:
        page.keyboard.press("Escape")
        page.evaluate("window.FastMDXDashboard.navigate('overview')")
        if width <= 800:
            page.locator("#problems-phone").click()
        elif page.locator("#side-expand").is_visible():
            page.locator("#side-expand").click()
        page.locator("#problems-tab").click()
        page.get_by_role("button", name="Inspect retained force evidence").first.click()
        page.locator("#instability-canvas[data-ready='true']").wait_for()
        expect(page.locator("#instability-ranking button")).to_have_count(4)
        page.keyboard.press("Escape")
        if width <= 800: expect(page.locator("#problems-phone")).to_be_focused()
    finally: page.context.close()


def test_switching_studies_discards_a_late_force_response(browser, site, tmp_path):
    page = _open(browser, site)
    try:
        page.locator("#instability-canvas[data-ready='true']").wait_for()
        pending = []
        page.route("**/api/instability?*", lambda route: pending.append(route))
        page.locator("#instability-refresh").click()
        page.wait_for_function("document.querySelector('#instability-status').textContent.includes('Reading retained')")
        other = tmp_path / "other"; (other / "simulation").mkdir(parents=True)
        (other / "manifest.json").write_text(json.dumps({"system": "Other", "phases": [{"name": "simulation", "status": "completed"}]}))
        page.evaluate("""async folder => {const data=await (await fetch('/api/explore/switch', {method:'POST', headers:{'content-type':'application/json'}, body:JSON.stringify({folder})})).json(); window.FastMDXDashboard.studyOpened(data, 'overview');}""", str(other))
        expect(page.locator("#instability-dialog")).to_be_hidden()
        for route in pending:
            try: route.fulfill(json={"study": site["key"], "available": True, "pdb": "invalid"})
            except sync_api.Error: pass
        page.wait_for_function("window.FastMDXInstability && window.FastMDXInstability.engine === null")
        page.unroute("**/api/instability?*")
        page.reload(wait_until="domcontentloaded")
        page.locator("#instability-open").click()
        expect(page.locator("#instability-status")).to_contain_text("Older runs did not record")
        assert not site["calls"]
    finally: page.context.close()


def test_late_evidence_cannot_reopen_a_closed_dialog(browser, site):
    page = _open(browser, site)
    try:
        page.locator("#instability-canvas[data-ready='true']").wait_for()
        pending = []
        page.route("**/api/instability?*", lambda route: pending.append(route))
        page.locator("#instability-refresh").click()
        page.wait_for_function("document.querySelector('#instability-status').textContent.includes('Reading retained')")
        page.keyboard.press("Escape")
        for route in pending:
            try: route.fulfill(json={"study": site["key"], "available": True, "pdb": "invalid"})
            except sync_api.Error: pass
        expect(page.locator("#instability-dialog")).to_be_hidden()
        assert page.evaluate("window.FastMDXInstability.engine === null")
        assert not site["calls"]
    finally: page.context.close()


def test_a_delayed_first_atom_selection_cannot_replace_the_last_choice(browser, site):
    page = _open(browser, site)
    try:
        page.locator("#instability-canvas[data-ready='true']").wait_for()
        # Engine methods are already queued. Delay the underlying renderer,
        # not a queued method waiting recursively on its own queue.
        page.evaluate("""() => {const renderer=window.FastMDXInstability.engine.plugin.builders.structure.representation;
          const original=renderer.addRepresentation.bind(renderer);window.pickEntered=0;
          renderer.addRepresentation=async (...args) => {const n=++window.pickEntered;
            if(n===1)await new Promise(resolve => window.releaseFirstPick=resolve);
            if(n===2)await new Promise(resolve => window.releaseSecondPick=resolve);
            return original(...args);};}""")
        page.locator("#instability-ranking button").first.click()
        page.wait_for_function("window.pickEntered === 1")
        page.locator("#instability-ranking button").last.click()
        page.evaluate("window.releaseFirstPick()")
        page.wait_for_function("window.pickEntered === 2")
        # No obsolete first choice is announced while the last renders.
        expect(page.locator("#instability-atom-note")).to_have_text("")
        page.evaluate("window.releaseSecondPick()")
        expect(page.locator("#instability-atom-note")).to_contain_text("GLY CA")
        atom = page.evaluate("() => {const e=window.FastMDXInstability.engine; return e.atoms(e.scene.selected)[0];}")
        assert atom["atom"] == "CA" and atom["icode"] == "B"
    finally: page.context.close()
