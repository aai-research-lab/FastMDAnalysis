"""⌘K grown into quick open, and every key on one sheet.

⌘K (Ctrl+K) finds a page with its key, every study of the workspace, a file
of the open study and a setting; ⌘P (Ctrl+P) the files alone. `@TRP23`
selects that residue in the Viewer and `:40` shows that frame there. The
study's pages are on ⌘1 to ⌘5, and on ⌥1 to ⌥5 (Alt) where a browser keeps
⌘1 for its tabs; ⌘B folds the sidebar and ⌘⌥B the side panel; ? lists every
key, the Viewer's with them.
"""

from __future__ import annotations

from pathlib import Path

import pytest

pytest.importorskip("playwright.sync_api")

from tests.test_the_cartoon_is_dssp_of_each_frame import _helical_study  # noqa: E402
from tests.test_the_workspace_says_its_studies import _study  # noqa: E402


def _serve(root: Path, look):
    from playwright.sync_api import sync_playwright

    from fastmdxplora.gui.server import start_dashboard_session

    session = start_dashboard_session(output=str(root), host="127.0.0.1", port=0)
    try:
        with sync_playwright() as pw:
            browser = pw.chromium.launch(args=["--enable-unsafe-swiftshader"])
            page = browser.new_page(viewport={"width": 1440, "height": 900})
            page.set_default_timeout(90000)
            errors: list[str] = []
            page.on("pageerror", lambda error: errors.append(str(error)))
            said = look(page, session.url)
            browser.close()
    finally:
        session.server.shutdown()
    assert errors == []
    return said


LINES = """() => [...document.querySelectorAll('#palette-list .palette-item')].map((li) => [
    li.closest('ul').querySelectorAll('.palette-group').length ? (() => {
        let g = li.previousElementSibling;
        while (g && !g.classList.contains('palette-group')) g = g.previousElementSibling;
        return g ? g.textContent : ''; })() : '',
    li.querySelector('.palette-label').textContent,
    [...li.querySelectorAll('.palette-keys kbd')].map((k) => k.textContent).join(' ')])"""


@pytest.fixture(scope="module")
def workspace(tmp_path_factory) -> Path:
    root = tmp_path_factory.mktemp("quick")
    _study(root / "ubiquitin", means={"rmsd": (0.1, 0.01)})
    _study(root / "lysozyme", system="1L63")
    return root


def test_k_finds_pages_with_their_keys_studies_files_and_settings(workspace):
    def look(page, url):
        page.goto(url + "#overview", wait_until="domcontentloaded")
        page.wait_for_function("() => window.FastMDXPalette && window.FastMDXKeys")
        page.wait_for_timeout(800)
        page.keyboard.press("Control+k")
        page.wait_for_selector("#palette:not([hidden])")
        page.wait_for_function(
            "() => [...document.querySelectorAll('#palette-list .palette-item')]"
            ".some((li) => li.textContent.includes('lysozyme'))")
        everything = page.evaluate(LINES)
        page.fill("#palette-input", "manifest")
        page.wait_for_timeout(200)
        files = page.evaluate(LINES)
        page.fill("#palette-input", "representation")
        page.wait_for_timeout(200)
        setting = page.evaluate(LINES)
        page.keyboard.press("Escape")
        page.keyboard.press("Control+p")
        page.wait_for_selector("#palette:not([hidden])")
        page.wait_for_timeout(500)
        only = page.evaluate(LINES)
        placeholder = page.get_attribute("#palette-input", "placeholder")
        page.keyboard.press("Escape")
        return {"everything": everything, "files": files, "setting": setting, "only": only,
                "placeholder": placeholder}

    said = _serve(workspace / "ubiquitin", look)
    pages = {line[1]: line[2] for line in said["everything"] if line[0] == "Page"}
    # In a browser's tab: the keys a tab leaves to the page.
    assert pages["Overview"] == "Alt 1" and pages["Files"] == "Alt 5"
    studies = [line[1] for line in said["everything"] if line[0] == "Study"]
    assert any("lysozyme" in s for s in studies) and any("ubiquitin" in s for s in studies)
    assert any(line[0] == "File" and line[1] == "manifest.json" for line in said["files"])
    assert any(line[0] == "Setting" and line[1] == "Representation" for line in said["setting"])
    assert said["only"] and {line[0] for line in said["only"]} == {"File"}
    assert said["placeholder"].startswith("A file of the open study")


def test_the_keys_open_pages_fold_the_columns_and_list_every_key(workspace):
    def look(page, url):
        page.goto(url + "#overview", wait_until="domcontentloaded")
        page.wait_for_function("() => window.FastMDXKeys && window.FastMDXFrame")
        page.wait_for_timeout(800)
        seen = {}
        for key, page_name in (("Alt+3", "analysis"), ("Alt+4", "report"), ("Alt+1", "overview")):
            page.keyboard.press(key)
            page.wait_for_function(f"() => document.documentElement.dataset.page === '{page_name}'")
            seen[key] = page_name
        folded = page.evaluate("() => document.body.classList.contains('sidebar-collapsed')")
        page.keyboard.press("Control+b")
        page.wait_for_timeout(300)
        after = page.evaluate("() => document.body.classList.contains('sidebar-collapsed')")
        page.keyboard.press("Control+b")
        page.wait_for_timeout(300)
        back = page.evaluate("() => document.body.classList.contains('sidebar-collapsed')")
        panel = page.evaluate("() => document.body.classList.contains('panel-collapsed')")
        page.keyboard.press("Control+Alt+b")
        page.wait_for_timeout(300)
        panel_after = page.evaluate("() => document.body.classList.contains('panel-collapsed')")
        page.keyboard.press("?")
        page.wait_for_selector("#viewer-help:not([hidden])")
        sheet = page.eval_on_selector_all(
            "#every-key > div", "rows => rows.map((r) => [r.querySelector('dt').textContent.trim(),"
            " r.querySelector('dd').textContent.trim()])")
        viewer = page.text_content("#viewer-keys")
        return {"seen": seen, "sidebar": [folded, after, back], "panel": [panel, panel_after],
                "sheet": sheet, "viewer": viewer}

    said = _serve(workspace / "ubiquitin", look)
    assert said["sidebar"][1] != said["sidebar"][0] and said["sidebar"][2] == said["sidebar"][0]
    assert said["panel"][1] != said["panel"][0]
    keys = dict(said["sheet"])
    assert "Go to anything" in keys["Ctrl K"] and "file" in keys["Ctrl P"]
    assert keys["Ctrl 3"] == "The study's Analysis (also Alt 3)"
    assert "Fold or show the side panel" == keys["Ctrl Alt B"]
    assert "Space" in said["viewer"]


def test_at_selects_a_residue_and_colon_shows_a_frame(tmp_path):
    study = _helical_study(tmp_path / "study", analyse=False)

    def look(page, url):
        page.goto(url + "#viewer", wait_until="domcontentloaded")
        page.wait_for_function(
            "() => window.FastMDXMoleculeViewer && window.FastMDXMoleculeViewer.STATE.engine"
            " && window.FastMDXMoleculeViewer.STATE.engine.structure()")
        page.wait_for_timeout(1000)
        page.keyboard.press("Control+k")
        page.fill("#palette-input", "@B:HIS63")
        page.wait_for_timeout(200)
        asked = page.evaluate(LINES)
        page.keyboard.press("Enter")
        page.wait_for_function(
            "() => (window.FastMDXMoleculeViewer.STATE.selection || {}).residues")
        chosen = page.evaluate("() => window.FastMDXMoleculeViewer.STATE.selection.residues")
        page.keyboard.press("Control+k")
        page.fill("#palette-input", "@LEU")
        page.wait_for_timeout(200)
        leucines = page.evaluate(LINES)
        page.fill("#palette-input", ":3")
        page.keyboard.press("Enter")
        page.wait_for_function("() => document.getElementById('traj-slider')"
                               " && document.getElementById('traj-slider').value === '3'")
        return {"asked": asked, "chosen": chosen, "leucines": leucines}

    said = _serve(study, look)
    assert said["asked"][0][1] == "Select B:HIS63 in the Viewer"
    assert said["chosen"] == [{"chain": "B", "resi": 63, "icode": "", "resn": "HIS"}]
    assert said["leucines"] and all(line[0] == "Residue" and "LEU" in line[1]
                                    for line in said["leucines"])
