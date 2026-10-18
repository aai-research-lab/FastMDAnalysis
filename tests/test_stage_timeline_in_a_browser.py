"""The stage molecule stays a saved snapshot while the study keeps running."""
from __future__ import annotations

import json
import shutil
from pathlib import Path

import pytest

sync_api = pytest.importorskip("playwright.sync_api")
expect = sync_api.expect
DATA = Path(__file__).parent / "data" / "stage_states"


@pytest.fixture(scope="module")
def browser():
    with sync_api.sync_playwright() as p:
        browser = p.chromium.launch(args=["--enable-unsafe-swiftshader"])
        yield browser
        browser.close()


@pytest.fixture
def site(tmp_path, monkeypatch):
    from fastmdxplora.gui.server import start_dashboard_session
    monkeypatch.setenv("FASTMDXPLORA_CONFIG_DIR", str(tmp_path / "config"))
    root=tmp_path/"study"; (root/"setup").mkdir(parents=True); (root/"simulation").mkdir()
    (root/"manifest.json").write_text(json.dumps({"system":"1UBQ","phases":[
        {"name":"setup","status":"completed"},{"name":"simulation","status":"completed"}]}))
    for name in ["input.pdb","prepared.pdb","solvated.pdb"]:
        shutil.copy(DATA/"topology.pdb",root/"setup"/name)
    shutil.copy(DATA/"topology.pdb",root/"simulation"/"topology.pdb")
    for name in ["state_minimized.xml","state_final.xml"]:
        shutil.copy(DATA/"state.xml",root/"simulation"/name)
    (root/"simulation"/"live_status.json").write_text(json.dumps({"status":"completed",
        "stage":"production","system_id":"1UBQ", "stage_states":{
            "minimization":"completed","nvt":"completed","npt":"completed","production":"completed"}}))
    session=start_dashboard_session(output=str(root),exploration_root=tmp_path,host="127.0.0.1",port=0)
    yield {"root":root,"session":session}
    session.stop()


def _open(browser,site,width=1440,theme="dark"):
    context=browser.new_context(viewport={"width":width,"height":900})
    context.add_init_script(f"localStorage.setItem('fmx.theme','{theme}')")
    page=context.new_page(); page.set_default_timeout(20000)
    page.errors=[]; page.on("pageerror",lambda e:page.errors.append(str(e)))
    page.goto(site["session"].url+"#overview",wait_until="domcontentloaded")
    page.locator('#stage-timeline [data-stage="minimized"]').wait_for()
    return page


@pytest.mark.parametrize("nested",[False,True])
def test_a_minimized_structure_is_the_states_coordinates_in_a_frozen_canvas(browser,site,monkeypatch,nested):
    if nested:
        child=site["root"]/"runs"/"replica";child.mkdir(parents=True)
        for folder in ["setup","simulation"]:shutil.copytree(site["root"]/folder,child/folder)
        shutil.copy(site["root"]/"manifest.json",child/"manifest.json")
        monkeypatch.setattr("fastmdxplora.gui.runs_together.run_shown",lambda root:child)
    page=_open(browser,site)
    try:
        page.click('#stage-timeline [data-stage="minimized"]')
        page.locator('#stage-structure-canvas[data-ready="true"]').wait_for()
        expect(page.locator("#stage-structure-title")).to_have_text("Minimized: retained structure")
        expect(page.locator("#stage-structure-source")).to_contain_text("state_minimized.xml")
        if nested:
            expect(page.locator("#stage-structure-source")).to_contain_text("runs/replica/simulation/state_minimized.xml")
            expect(page.locator("#stage-timeline-note")).to_contain_text("For runs/replica.")
        point=page.evaluate("""() => {
            const engine=window.FastMDXStageTimeline.engine;
            const structure=engine.structure(); const unit=structure.units[0];
            const point=[0,0,0];unit.conformation.position(unit.elements[0],point);return point;
        }""")
        assert point == pytest.approx([3,4,5],abs=.001)
        page.evaluate("""() => window.FastMDXDashboard.applyStatus({
            study:window.FastMDXDashboard.state.studyShown,
            status:{stage:'production',status:'running',system_id:'1UBQ'}})""")
        expect(page.locator("#stage-structure-title")).to_have_text("Minimized: retained structure")
        expect(page.locator("#stage-structure-canvas")).to_have_attribute("data-ready","true")
        page.get_by_role("button",name="Return to study",exact=True).click()
        expect(page.locator("#stage-structure-dialog")).to_be_hidden()
        page.wait_for_function("window.FastMDXStageTimeline.engine === null")
        expect(page.locator('#stage-timeline [data-stage="minimized"]')).to_be_focused()
        assert page.errors == []
    finally: page.context.close()


def test_a_completed_stage_without_coordinates_explains_the_missing_snapshot(browser,site):
    page=_open(browser,site,390)
    try:
        page.click('[data-stage="equilibrated"]')
        expect(page.locator("#stage-structure-note")).to_contain_text("not retained")
        assert page.locator("#stage-structure-canvas canvas").count()==0
        expect(page.locator("#stage-structure-solvent")).to_be_disabled()
        assert page.locator(".stage-structure-box").bounding_box()["height"] < 450
        assert page.get_by_role("button",name="Return to study",exact=True).bounding_box()["width"] >= 100
        page.keyboard.press("Escape")
        expect(page.locator("#stage-structure-dialog")).to_be_hidden()
        assert page.errors == []
    finally: page.context.close()


@pytest.mark.parametrize("kind",["ligand","water"])
def test_retained_atoms_are_rendered_without_current_study_classification(browser,site,kind):
    target=site["root"]/"setup"/"input.pdb"
    text=target.read_text()
    if kind=="ligand":
        text=text.replace("ATOM  ","HETATM").replace("GLY","LIG")
    else:
        text=text.replace("END\n","HETATM    5  O   HOH W   2       8.000   8.000   8.000  1.00  0.00           O\nEND\n")
    target.write_text(text)
    page=_open(browser,site)
    try:
        page.click('[data-stage="input"]')
        page.locator('#stage-structure-canvas[data-ready="true"]').wait_for()
        if kind=="water":
            page.check("#stage-structure-solvent")
            expect(page.locator("#stage-structure-note")).to_contain_text("5 atoms")
        count=page.evaluate("window.FastMDXStageTimeline.engine.components.retained.cell.obj.data.elementCount")
        assert count==(5 if kind=="water" else 4)
        assert page.errors==[]
    finally: page.context.close()


def test_an_invalid_saved_state_is_not_replaced_by_the_current_molecule(browser,site):
    (site["root"]/"simulation"/"state_minimized.xml").write_text('<State type="State" version="1"/>')
    page=_open(browser,site)
    try:
        page.click('[data-stage="minimized"]')
        expect(page.locator("#stage-structure-note")).to_contain_text("no readable positions")
        assert page.locator("#stage-structure-canvas canvas").count()==0
        assert page.errors == []
    finally: page.context.close()


def test_closing_a_slow_snapshot_prevents_a_late_render(browser,site):
    page=_open(browser,site)
    try:
        held=[]
        page.route("**/structure/stage.pdb?*",lambda r:held.append(r))
        with page.expect_request("**/structure/stage.pdb?*"):
            page.click('[data-stage="minimized"]')
        page.keyboard.press("Escape")
        expect(page.locator("#stage-structure-dialog")).to_be_hidden()
        held[0].fulfill(content_type="chemical/x-pdb",body=(DATA/"topology.pdb").read_text())
        page.evaluate("() => new Promise(requestAnimationFrame)")
        assert page.locator("#stage-structure-canvas canvas").count()==0
        assert page.evaluate("window.FastMDXStageTimeline.engine === null")
        assert page.errors == []
    finally: page.context.close()


def test_snapshot_dialog_is_not_restored_on_reload(browser,site):
    page=_open(browser,site)
    try:
        page.click('[data-stage="minimized"]')
        expect(page.locator("#stage-structure-dialog")).to_be_visible()
        page.reload(wait_until="domcontentloaded")
        page.locator('#stage-timeline [data-stage="minimized"]').wait_for()
        expect(page.locator("#stage-structure-dialog")).to_be_hidden()
        assert page.errors == []
    finally: page.context.close()


def test_a_late_snapshot_cannot_reopen_after_switching_studies(browser,site):
    other=site["root"].parent/"second"
    shutil.copytree(site["root"],other)
    (other/"simulation"/"state_minimized.xml").unlink()
    page=_open(browser,site)
    try:
        held=[]
        page.route("**/structure/stage.pdb?*",lambda route:held.append(route))
        with page.expect_request("**/structure/stage.pdb?*"):
            page.click('[data-stage="minimized"]')
        page.evaluate("""async path => {
            const response=await fetch('/api/explore/switch', {method:'POST',
                headers:{'content-type':'application/json'},body:JSON.stringify({folder:path})});
            const answer=await response.json();
            if (!answer.ok) throw new Error(JSON.stringify(answer));
            window.FastMDXDashboard.studyOpened(answer,'overview');
        }""",str(other))
        expect(page.locator("#stage-structure-dialog")).to_be_hidden()
        expect(page.locator('[data-stage="minimized"]')).to_have_class("stage-unavailable")
        held[0].fulfill(content_type="chemical/x-pdb",body=(DATA/"topology.pdb").read_text())
        page.evaluate("() => new Promise(requestAnimationFrame)")
        assert page.locator("#stage-structure-canvas canvas").count()==0
        assert page.evaluate("window.FastMDXStageTimeline.engine === null")
        assert page.errors==[]
    finally: page.context.close()


def test_closing_returns_focus_after_the_timeline_updates(browser,site):
    page=_open(browser,site)
    try:
        page.click('[data-stage="minimized"]')
        page.locator('#stage-structure-canvas[data-ready="true"]').wait_for()
        (site["root"]/"simulation"/"live_status.json").write_text(json.dumps({"stage":"minimization",
            "status":"running","stage_states":{"minimization":"current"}}))
        page.evaluate("window.FastMDXStageTimeline.refresh()")
        expect(page.locator('[data-stage="minimized"]')).to_have_attribute("data-state","current")
        page.get_by_role("button",name="Return to study",exact=True).click()
        expect(page.locator('[data-stage="minimized"]')).to_be_focused()
        assert page.errors==[]
    finally: page.context.close()


@pytest.mark.parametrize("width",[390,980,1600])
@pytest.mark.parametrize("theme",["light","dark"])
def test_timeline_and_snapshot_fit_and_close_by_keyboard(browser,site,width,theme):
    page=_open(browser,site,width,theme)
    try:
        row=page.locator('.stage-timeline-row').bounding_box()
        timeline_note=page.locator('#stage-timeline-note').bounding_box()
        assert timeline_note['y'] >= row['y']+row['height']-1
        assert row['width'] >= width-260
        stage=page.locator('[data-stage="minimized"]'); stage.focus(); page.keyboard.press("Enter")
        page.locator('#stage-structure-canvas[data-ready="true"]').wait_for()
        assert page.evaluate("document.documentElement.scrollWidth <= innerWidth")
        bounds=page.locator(".stage-structure-box").bounding_box()
        assert bounds["x"]>=0 and bounds["x"]+bounds["width"]<=width+1
        assert bounds["y"]>=0 and bounds["y"]+bounds["height"]<=901
        expect(page.get_by_role("button",name="Return to study",exact=True)).to_be_in_viewport()
        button=page.get_by_role("button",name="Return to study",exact=True).bounding_box()
        assert button["width"] >= 100 and button["height"] < 50
        page.keyboard.press("Escape")
        expect(stage).to_be_focused()
        assert page.errors == []
    finally: page.context.close()
