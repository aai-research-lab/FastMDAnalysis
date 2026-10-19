"""The Analysis page as cells, its plots each zoomed, moved and saved.

Each analysis was a heading over its figures; what it determined was in the
table at the top, what it was run with behind a chip, and Analyze again ran
every analysis. Each is now a cell: under its heading, its values as the
table says them, what it was run with (the frames, the atoms, its options)
and Analyze again with it alone ticked. Every plot zooms by a drag across
it, moves with Shift held, shows all again by a double click or its
button, and is saved as SVG, PNG or CSV of what it shows. A time range
dragged across a plot on the simulation's clock is every such plot's, said
at the top of the page and marked on the Viewer's timeline. A map
(secondary structure over time, a matrix, a landscape) is drawn as its
cells, not as a picture of them, and a run under a nanosecond is read in
picoseconds on every plot.
"""

from __future__ import annotations

from pathlib import Path

import pytest


@pytest.fixture(scope="module")
def study(tmp_path_factory) -> Path:
    from tests.test_the_cartoon_is_dssp_of_each_frame import _helical_study

    root = _helical_study(tmp_path_factory.mktemp("cells") / "study", analyse=False)
    from fastmdxplora.analysis import AnalysisOrchestrator

    AnalysisOrchestrator(str(root / "simulation" / "production.dcd"),
                         str(root / "simulation" / "trajectory_topology.pdb"),
                         output_dir=str(root / "analysis")).run(include=["rmsd", "rg", "ss"])
    return root


LOOK = """() => {
    const heads = Object.fromEntries([...document.querySelectorAll('.analysis-cell-head')]
        .map((h) => [h.dataset.analysis, h.textContent]));
    const chart = (sel) => document.querySelector(sel);
    const ss = chart('.figure-chart[data-figure-kind="states"]');
    return {
        heads,
        rg: (chart('.series-chart[data-analysis="rg"]') || {}).dataset?.view,
        rmsd: (chart('.series-chart[data-analysis="rmsd"]') || {}).dataset?.view,
        ss: ss ? ss.dataset.view : null,
        ssCells: ss ? ss.querySelector('svg').dataset.cells : null,
        ssImages: ss ? ss.querySelectorAll('svg image').length : null,
        ssWords: ss ? [...ss.querySelectorAll('svg text')].map((t) => t.textContent) : [],
        range: document.getElementById('analysis-cell-range')
            && !document.getElementById('analysis-cell-range').hidden
            ? document.getElementById('analysis-cell-range').textContent : '',
    };
}"""


def test_each_analysis_is_a_cell_whose_plots_zoom_share_a_range_and_save(study, tmp_path):
    pytest.importorskip("playwright.sync_api")
    from playwright.sync_api import sync_playwright

    from fastmdxplora.gui.server import start_dashboard_session

    session = start_dashboard_session(output=str(study), host="127.0.0.1", port=0)
    try:
        with sync_playwright() as pw:
            browser = pw.chromium.launch(args=["--enable-unsafe-swiftshader"])
            page = browser.new_page(viewport={"width": 1440, "height": 900}, accept_downloads=True)
            page.set_default_timeout(60000)
            errors: list[str] = []
            page.on("pageerror", lambda error: errors.append(str(error)))
            page.goto(session.url + "#analysis", wait_until="domcontentloaded")
            page.wait_for_selector('.analysis-cell-head[data-analysis="rmsd"] .cell-values')
            page.wait_for_selector('.series-chart[data-analysis="rmsd"] .plot-tools')
            page.wait_for_selector('.figure-chart[data-figure-kind="states"] .plot-tools')
            before = page.evaluate(LOOK)

            # Analyze again from the RMSD's cell: it alone ticked.
            page.click('.analysis-cell-head[data-analysis="rmsd"] .cell-again')
            page.wait_for_selector("#analysis-again-ask .again-choices")
            ticked = page.evaluate("""() => [...document.querySelectorAll(
                '#analysis-again-ask .again-choices input:checked')].map((i) => i.value)""")
            page.click("#analysis-again-ask .fix-cancel")

            # A range dragged across the RMSD: every plot on the clock shows it.
            svg = page.locator('.series-chart[data-analysis="rmsd"] svg')
            svg.scroll_into_view_if_needed()
            box = svg.bounding_box()
            page.mouse.move(box["x"] + box["width"] * 0.35, box["y"] + 100)
            page.mouse.down()
            page.mouse.move(box["x"] + box["width"] * 0.75, box["y"] + 110, steps=6)
            page.mouse.up()
            page.wait_for_function("() => document.querySelector("
                                   "'.series-chart[data-analysis=\"rg\"]').dataset.view === 'range'")
            ranged = page.evaluate(LOOK)
            with page.expect_download() as got:
                page.click('.series-chart[data-analysis="rmsd"] .plot-tool:text("CSV")')
            csv = Path(got.value.path()).read_text(encoding="utf-8").splitlines()
            with page.expect_download() as got:
                page.click('.figure-chart[data-figure-kind="states"] .plot-tool:text("SVG")')
            svg_file = (got.value.suggested_filename, Path(got.value.path()).read_text(encoding="utf-8"))
            with page.expect_download() as got:
                page.click('.series-chart[data-analysis="rg"] .plot-tool:text("PNG")')
            png = (got.value.suggested_filename, Path(got.value.path()).read_bytes()[:8])
            # A map of more runs of a colour than it draws as cells: an image of
            # them, and said so.
            many = page.evaluate("""() => {
                const host = document.createElement('div');
                host.style.width = '600px';
                document.querySelector('.page[data-page="analysis"]').appendChild(host);
                const frames = 400, residues = 300;
                const x = [...Array(frames).keys()].map((i) => 0.001 * (i + 1));
                window.FastMDXFigures.plot(host, {
                    ok: true, kind: 'states', title: 'Secondary structure', x: x,
                    frames: [...Array(frames).keys()], x_label: 'Time (ns)', y_label: 'Residue',
                    residues: [...Array(residues).keys()].map((i) => String(i + 1)),
                    states: x.map((_, i) => [...Array(residues).keys()]
                        .map((r) => ((i + r) % 2 ? 'H' : 'E')).join('')),
                    classes: [{code: 'H', label: 'Helix'}, {code: 'E', label: 'Strand'}],
                    linked: false});
                const image = host.querySelector('svg image');
                const said = [host.querySelector('svg').dataset.cells,
                    [...host.querySelectorAll('svg text')].map((t) => t.textContent),
                    image.querySelector('title').textContent, getComputedStyle(image).imageRendering];
                host.remove();
                return said;
            }""")
            # Show all, at the top of the page.
            page.click("#analysis-cell-range .plot-tool")
            page.wait_for_function("() => document.querySelector("
                                   "'.series-chart[data-analysis=\"rg\"]').dataset.view === ''")
            after = page.evaluate(LOOK)
            browser.close()
    finally:
        session.server.shutdown()
    assert errors == []
    rmsd = before["heads"]["rmsd"]
    assert rmsd.startswith("Mean") and "6 frames" in rmsd and "Analyze again" in rmsd
    assert "Helix fraction" in before["heads"]["ss"]
    assert ticked == ["rmsd"]
    # A map drawn as its cells.
    assert before["ssCells"].isdigit() and int(before["ssCells"]) > 0 and before["ssImages"] == 0
    # This run kept no clock: its plots are on frame numbers.
    assert "Frame" in before["ssWords"]
    assert (before["rg"], before["ss"]) == ("", "")
    assert (ranged["rmsd"], ranged["rg"], ranged["ss"]) == ("range", "range", "range")
    assert ranged["range"].startswith("Every plot on the simulation's clock shows ")
    assert csv[0] == "frame,rmsd (nm)"
    shown = [float(line.split(",")[0]) for line in csv[1:]]
    assert 0 < len(shown) < 6, csv
    assert svg_file[0].endswith(".svg") and svg_file[1].startswith("<svg") and "xmlns" in svg_file[1]
    assert png[0].endswith(".png") and png[1] == b"\x89PNG\r\n\x1a\n"
    assert (after["rmsd"], after["rg"], after["ss"], after["range"]) == ("", "", "", "")
    # A run under a nanosecond (0.4 ns here) read in picoseconds.
    assert many[0] == "image" and "Time (ps)" in many[1] and "400" in many[1]
    assert many[2].startswith("Drawn as a picture of its 400 by 300 cells") and many[3] == "pixelated"


def _page(pw, session):
    browser = pw.chromium.launch(args=["--enable-unsafe-swiftshader"])
    page = browser.new_page(viewport={"width": 1440, "height": 900}, accept_downloads=True)
    page.set_default_timeout(60000)
    errors: list[str] = []
    page.on("pageerror", lambda error: errors.append(str(error)))
    return browser, page, errors


PROJECTION = """() => {
    const host = document.createElement('div');
    host.id = 'a-projection';
    host.style.width = '600px';
    document.querySelector('.page[data-page="analysis"]').appendChild(host);
    const n = 100;
    window.FastMDXFigures.plot(host, {
        ok: true, kind: 'projection', title: 'Projection', axes: ['PC1', 'PC2'], unit: '',
        x: [...Array(n).keys()].map((i) => i % 10), y: [...Array(n).keys()].map((i) => Math.floor(i / 10)),
        t: [...Array(n).keys()], frames: [...Array(n).keys()], t_label: 'Frame', linked: false});
    host.scrollIntoView();
}"""


def _csv(page, selector):
    with page.expect_download() as got:
        page.click(selector + ' .plot-tool:text("CSV")')
    lines = Path(got.value.path()).read_text(encoding="utf-8").splitlines()
    return lines[0], [[float(cell) for cell in line.split(",")[2:]] for line in lines[1:]]


def _drag(page, selector, a, b, shift=False):
    box = page.locator(selector + " svg").bounding_box()
    if shift:
        page.keyboard.down("Shift")
    page.mouse.move(box["x"] + box["width"] * a[0], box["y"] + box["height"] * a[1])
    page.mouse.down()
    page.mouse.move(box["x"] + box["width"] * b[0], box["y"] + box["height"] * b[1], steps=6)
    page.mouse.up()
    if shift:
        page.keyboard.up("Shift")


def test_a_map_zooms_on_both_axes_moves_and_shows_all_again(study):
    pytest.importorskip("playwright.sync_api")
    from playwright.sync_api import sync_playwright

    from fastmdxplora.gui.server import start_dashboard_session

    session = start_dashboard_session(output=str(study), host="127.0.0.1", port=0)
    try:
        with sync_playwright() as pw:
            browser, page, errors = _page(pw, session)
            page.goto(session.url + "#analysis", wait_until="domcontentloaded")
            page.wait_for_selector('.analysis-cell-head[data-analysis="rmsd"]')
            page.evaluate(PROJECTION)
            host = "#a-projection"
            head, every = _csv(page, host)
            _drag(page, host, (0.3, 0.3), (0.6, 0.6))
            page.wait_for_function("() => document.querySelector('#a-projection').dataset.view === 'zoomed'")
            _, zoomed = _csv(page, host)
            _drag(page, host, (0.6, 0.5), (0.3, 0.5), shift=True)
            page.wait_for_timeout(300)
            moved_view = page.evaluate("() => document.querySelector('#a-projection').dataset.view")
            _, moved = _csv(page, host)
            box = page.locator(host + " svg").bounding_box()
            page.mouse.dblclick(box["x"] + box["width"] * 0.45, box["y"] + box["height"] * 0.5)
            page.wait_for_function("() => document.querySelector('#a-projection').dataset.view === ''")
            _, again = _csv(page, host)
            browser.close()
    finally:
        session.server.shutdown()
    assert errors == []
    assert head == "frame,frame_shown,PC1,PC2" and len(every) == 100 and len(again) == 100
    # Less of it on both axes: neither every x nor every y.
    assert 0 < len(zoomed) < 100
    assert len({x for x, _ in zoomed}) < 10 and len({y for _, y in zoomed}) < 10
    # Moved to the right of what was zoomed, and no higher or lower.
    assert moved_view == "zoomed"
    assert min(x for x, _ in moved) > min(x for x, _ in zoomed)
    assert {y for _, y in moved} == {y for _, y in zoomed}


def test_the_range_is_marked_on_the_viewer_s_timeline(study):
    pytest.importorskip("playwright.sync_api")
    from playwright.sync_api import sync_playwright

    from fastmdxplora.gui.server import start_dashboard_session

    session = start_dashboard_session(output=str(study), host="127.0.0.1", port=0)
    try:
        with sync_playwright() as pw:
            browser, page, errors = _page(pw, session)
            page.goto(session.url + "#analysis", wait_until="domcontentloaded")
            page.wait_for_selector('.series-chart[data-analysis="rmsd"] .plot-tools')
            _drag(page, '.series-chart[data-analysis="rmsd"]', (0.4, 0.4), (0.8, 0.45))
            page.wait_for_selector("#analysis-cell-range:not([hidden])")
            page.click('[data-view-link="viewer"]')
            page.wait_for_selector("#timeline-ruler svg .timeline-range")
            marked = page.evaluate("""() => [document.getElementById('frame-series').dataset.range,
                parseFloat(getComputedStyle(document.getElementById('frame-series'))
                    .getPropertyValue('--range-width'))]""")
            page.evaluate("() => window.FastMDXPlotTools.clearRange()")
            page.wait_for_function("() => !document.querySelector('#timeline-ruler svg .timeline-range')")
            cleared = page.evaluate("() => document.getElementById('frame-series').dataset.range || ''")
            browser.close()
    finally:
        session.server.shutdown()
    assert errors == []
    assert marked[0] == "on" and marked[1] > 1
    assert cleared == ""


@pytest.mark.parametrize("public", [False, True], ids=["loopback", "public"])
def test_beyond_loopback_no_cell_offers_analyze_again(study, public):
    """Analyze again POSTs, and a GUI served beyond this computer refuses
    it; each cell offered it there, only to be refused. The page offers it
    where it is answered."""
    pytest.importorskip("playwright.sync_api")
    from playwright.sync_api import sync_playwright

    from fastmdxplora.gui.exploration import DashboardRuntime
    from tests.test_a_public_dashboard_answers_only_what_it_lists import _serving

    runtime = DashboardRuntime(workspace_root=study, exploration_root=study.parent,
                               active_root=study)
    with _serving(runtime, allow_control=not public) as url, sync_playwright() as pw:
        browser, page, errors = _page(pw, None)
        page.goto(url + "#analysis", wait_until="domcontentloaded")
        page.wait_for_selector('.analysis-cell-head[data-analysis="rmsd"] .cell-values')
        page.wait_for_function("() => window.FastMDXDashboard && "
                               "window.FastMDXDashboard.state.appState.active_run")
        page.wait_for_timeout(500)
        offered = page.evaluate("""() => ({
            bound: document.body.dataset.beyondLoopback,
            cells: [...document.querySelectorAll('.cell-again')].filter((b) => !b.hidden).length,
            page: !document.getElementById('analysis-again').hidden,
        })""")
        browser.close()
    assert errors == []
    assert offered["bound"] == ("true" if public else "false")
    if public:
        assert offered == {"bound": "true", "cells": 0, "page": False}
    else:
        assert offered["cells"] >= 3 and offered["page"] is True
