"""A record the study pages did not read is said not read, not not determined.

The Overview and the Analysis page read only the study's own records: one
linked out of the study is refused (`record_access`). The mean whose record
was refused was badged "Not determined" and counted with the means not
determined ("0 of 1 mean determined"), though whether it was determined is
not known. And nothing on either page named a refused record: an analysis
manifest linked out quietly turned the analyses' times into frame numbers.

A refused mean is now Not read, counted apart, and each page names the
records it did not read, with what a clock not read does to the times.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from tests.test_the_study_pages_read_only_the_study import (
    CLOCKS,
    OUTSIDE,
    _ask,
    _link,
    _study,
    _write,
)


def _linked(tmp_path: Path, relative: str) -> Path:
    """The study with ``relative`` moved outside it and linked back."""
    root = _study(tmp_path)
    if relative == "analysis/rmsd/options.json":
        content = OUTSIDE["determined"]
    else:
        content = next(content for path, content in CLOCKS.values() if path == relative)
        if (root / relative).exists():
            (root / relative).unlink()
    outside = _write(tmp_path, "elsewhere/" + Path(relative).name, content)
    _link(root / relative, outside)
    return root


@pytest.mark.parametrize("public", [False, True], ids=["loopback", "public"])
def test_a_mean_not_read_is_counted_apart(tmp_path, public):
    root = _linked(tmp_path, "analysis/rmsd/options.json")

    overview = json.loads(_ask(root, tmp_path, "/api/overview", public=public)[1])
    analyses = json.loads(_ask(root, tmp_path, "/api/analysis-overview", public=public)[1])
    results = json.loads(_ask(root, tmp_path, "/api/results", public=public)[1])

    tiles = {tile["analysis"]: tile for tile in overview["tiles"]}
    assert tiles["rmsd"]["unavailable"] is True and tiles["rmsd"]["determined"] is None
    # Rg determined, of the one mean read; RMSD not read, not counted as a
    # mean not determined.
    assert (overview["means"], overview["determined"], overview["not_read"]) == (1, 1, 1)
    rows = {row["analysis"]: row for row in analyses["rows"]}
    assert rows["rmsd"]["quantities"][0]["unavailable"] is True
    metrics = {row["metric"]: row for row in results["metric_rows"]}
    assert metrics["RMSD"]["status"] == "Not read"
    assert metrics["Radius of gyration"]["status"] == "Determined"


# ---------------------------------------------------------------------------
# As the pages show it.
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def browser():
    pytest.importorskip("playwright.sync_api")
    from playwright.sync_api import sync_playwright

    with sync_playwright() as pw:
        launched = pw.chromium.launch(args=["--enable-unsafe-swiftshader"])
        yield launched
        launched.close()


def _pages(browser, root: Path) -> dict:
    """What the Overview and the Analysis page say of the study."""
    from fastmdxplora.gui.server import start_dashboard_session

    session = start_dashboard_session(output=str(root), host="127.0.0.1", port=0)
    page = browser.new_page(viewport={"width": 1440, "height": 900})
    page.set_default_timeout(60000)
    errors: list[str] = []
    page.on("pageerror", lambda error: errors.append(str(error)))
    try:
        page.goto(session.url + "#overview", wait_until="domcontentloaded")
        page.wait_for_selector("#overview-tiles .overview-tile")
        page.wait_for_selector("#overview-not-read")
        said = page.evaluate("""() => ({
            badges: Object.fromEntries([...document.querySelectorAll('#overview-tiles .overview-tile')]
                .map((t) => [t.dataset.tileAnalysis, t.querySelector('.tile-badge').textContent])),
            why: Object.fromEntries([...document.querySelectorAll('#overview-tiles .overview-tile')]
                .map((t) => [t.dataset.tileAnalysis, (t.querySelector('.tile-why') || {}).textContent || ''])),
            counted: document.getElementById('overview-results-said').textContent,
            line: document.getElementById('overview-not-read').textContent,
        })""")
        page.goto(session.url + "#analysis", wait_until="domcontentloaded")
        page.wait_for_function(
            "() => (document.getElementById('analysis-results-summary') || {}).textContent")
        said.update(page.evaluate("""() => ({
            summary: document.getElementById('analysis-results-summary').textContent,
            chips: Object.fromEntries([...document.querySelectorAll(
                '#analysis-results-table tr[data-analysis]')].map((row) => [
                row.dataset.analysis,
                (row.querySelector('.analysis-results-status') || {}).textContent || ''])),
        })"""))
    finally:
        page.close()
        session.stop()
    assert not errors, errors
    return said


def test_the_pages_say_a_mean_not_read(browser, tmp_path):
    said = _pages(browser, _linked(tmp_path, "analysis/rmsd/options.json"))

    assert said["badges"] == {"rmsd": "Not read", "rg": "Determined"}
    assert said["why"]["rmsd"].startswith("Not read. analysis/rmsd/options.json was not read")
    assert said["counted"].startswith("1 of 1 mean determined, 1 not read · ")
    assert said["line"].startswith("Not read: analysis/rmsd/options.json. It leads outside "
                                   "the study")
    assert said["chips"] == {"rmsd": "Not read", "rg": "Determined"}
    assert said["summary"].startswith("1 determined, 1 not read. ")
    assert "Not read: analysis/rmsd/options.json." in said["summary"]


@pytest.mark.parametrize("record", ["analysis/analysis_manifest.json",
                                    "simulation/simulation_parameters.json"])
def test_the_pages_say_a_clock_not_read(browser, tmp_path, record):
    said = _pages(browser, _linked(tmp_path, record))

    for text in (said["line"], said["summary"]):
        assert f"Not read: {record}." in text
        clock = "their frames are numbered, not timed" in text
        assert clock is (record == "analysis/analysis_manifest.json")
    # The mean itself was read.
    assert said["badges"] == {"rg": "Determined"} and said["chips"]["rg"] == "Determined"
