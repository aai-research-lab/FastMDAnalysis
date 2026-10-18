"""The first round of the Integrated MD Study Environment section, as seen.

A study of several runs said "Unknown" in the sidebar and on its Overview and
"completed" in All studies; its runs' dots were orange when completed and
their names cut to "random_see..."; its Analysis page said "no analyses yet"
above runs of fourteen each; the Overview's "Reset zoom" reset no zoom; a
chain numbered from 70 was drawn all red under Spectrum.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

pytest.importorskip("playwright.sync_api")

from tests.test_a_study_of_runs_is_shown_as_one import _a_sweep_study  # noqa: E402

DONE = {"stage": "completed", "current_step": 500, "total_planned_steps": 500}


def _serve(root: Path, look, *, width: int = 1440, colour_scheme: str = "dark"):
    from playwright.sync_api import sync_playwright

    from fastmdxplora.gui.server import start_dashboard_session

    session = start_dashboard_session(output=str(root), host="127.0.0.1", port=0)
    try:
        with sync_playwright() as pw:
            browser = pw.chromium.launch(args=["--enable-unsafe-swiftshader"])
            page = browser.new_page(viewport={"width": width, "height": 900},
                                    color_scheme=colour_scheme)
            page.set_default_timeout(60000)
            errors: list[str] = []
            page.on("pageerror", lambda error: errors.append(str(error)))
            said = look(page, session.url)
            browser.close()
    finally:
        session.server.shutdown()
    assert errors == []
    return said


@pytest.fixture(scope="module")
def finished_sweep(tmp_path_factory) -> Path:
    root = _a_sweep_study(tmp_path_factory.mktemp("ide"))
    for folder in sorted((root / "runs").iterdir()):
        (folder / "simulation" / "live_status.json").write_text(json.dumps(DONE), encoding="utf-8")
    return root


@pytest.mark.parametrize("width,scheme", [(1440, "dark"), (900, "light")])
def test_a_finished_study_of_runs_is_completed_everywhere(finished_sweep, width, scheme):
    def look(page, url):
        page.goto(url + "#overview", wait_until="domcontentloaded")
        page.wait_for_selector("#study-runs:not([hidden]) .study-run-row")
        page.wait_for_function(
            "() => document.getElementById('topbar-status-text').textContent !== 'Ready'")
        page.wait_for_timeout(500)
        said = {
            "word": page.text_content("#topbar-status-text"),
            "dot": page.get_attribute("#topbar-status-dot", "class"),
            "verdict": page.text_content("#verdict-word"),
            "runs": page.eval_on_selector_all(
                "#study-runs .study-run-row",
                "rs => rs.map(r => [r.querySelector('.study-running-name').textContent,"
                " r.querySelector('.study-running-name').title.split('\\n')[0],"
                " r.querySelector('.status-dot').className,"
                " [r.querySelector('.study-running-name').scrollWidth,"
                "  r.querySelector('.study-running-name').clientWidth]])"),
        }
        page.goto(url + "#studies", wait_until="domcontentloaded")
        page.wait_for_selector(".study-state, .studies-state")
        said["card"] = page.eval_on_selector_all(
            ".study-state", "s => s.map(x => x.textContent.trim())")
        return said

    said = _serve(finished_sweep, look, width=width, colour_scheme=scheme)
    assert said["word"] == "Completed" and said["verdict"] == "Completed"
    assert "status-dot-completed" in said["dot"]
    assert said["card"] == ["completed"]
    names = [row[0] for row in said["runs"]]
    assert names == ["temperature 300 K", "temperature 310 K"]
    # The folder's name on hover; each dot green; each name whole.
    assert all(row[1].startswith("s") and "temperature" in row[1] for row in said["runs"])
    assert all("status-dot-completed" in row[2] for row in said["runs"])
    assert all(row[3][0] <= row[3][1] for row in said["runs"]), said["runs"]
