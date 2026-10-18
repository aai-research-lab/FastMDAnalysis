"""The pages say what was not read, never that nothing was recorded.

`/api/metrics`, `/api/residue-values`, `/api/report` and `/api/fixes` name
a record of the study they did not read (a link out of it) and say why; the
pages ignored it. The charts were drawn empty under "The run recorded no
sample", the Viewer's "Coloured by" list lost its results without a word,
the Report page said "No report yet" above the line saying the report was
not read, and the remedy card gave the words of a Manifest read through a
link. Each now says the record was not read, and why.
"""

from __future__ import annotations

import json
import urllib.request
from pathlib import Path

import pytest

from tests.test_a_remedy_is_not_priced_from_a_record_not_read import _linked_out, _stopped
from tests.test_the_report_and_its_neighbours_read_only_the_study import _with_its_report

MARK = "OUTSIDE-PAGE-not-for-any-page"


@pytest.fixture(scope="module")
def browser():
    pytest.importorskip("playwright.sync_api")
    from playwright.sync_api import sync_playwright

    with sync_playwright() as pw:
        launched = pw.chromium.launch(args=["--enable-unsafe-swiftshader"])
        yield launched
        launched.close()


def _linked(tmp_path: Path) -> Path:
    """A study carried on whose report, per-residue result and piece's
    charts are links out of it."""
    root = _with_its_report(tmp_path)
    for record in ("report/report.md", "analysis/rmsf/options.json",
                   "segment-001/simulation/live_metrics.csv"):
        _linked_out(root, tmp_path, record)
    return root


def _said(browser, root: Path, page_name: str, ready: str, read: str) -> tuple[dict, dict]:
    """What the page ``page_name`` says once ``ready`` holds, and the
    answers its routes gave."""
    from fastmdxplora.gui.server import start_dashboard_session

    session = start_dashboard_session(output=str(root), host="127.0.0.1", port=0)
    page = browser.new_page(viewport={"width": 1440, "height": 900})
    page.set_default_timeout(60000)
    errors: list[str] = []
    page.on("pageerror", lambda error: errors.append(str(error)))
    try:
        page.goto(session.url + "#" + page_name, wait_until="domcontentloaded")
        page.wait_for_function(ready)
        said = page.evaluate(read)
        answers = {route: json.loads(urllib.request.urlopen(
            session.url.rstrip("/") + route, timeout=30).read())
            for route in ("/api/metrics", "/api/residue-values", "/api/report", "/api/fixes")}
    finally:
        page.close()
        session.stop()
    assert not errors, errors
    return said, answers


def test_the_charts_say_their_record_was_not_read(browser, tmp_path):
    said, answers = _said(
        browser, _linked(tmp_path), "overview",
        "() => window.FastMDXCharts && (window.FastMDXCharts.drawn('temperature') || {}).said",
        """() => ({empty: document.getElementById('chart-empty').textContent,
                   drawn: window.FastMDXCharts.drawn('temperature').said})""")
    not_read = answers["/api/metrics"]["not_read"]
    assert not_read.startswith("segment-001/simulation/live_metrics.csv was not read")
    assert said == {"empty": not_read, "drawn": "not read"}


def test_the_viewer_says_per_residue_results_were_not_read(browser, tmp_path):
    said, answers = _said(
        browser, _linked(tmp_path), "viewer",
        "() => !document.getElementById('viewer-colour-said').hidden",
        """() => ({note: document.getElementById('viewer-colour-said').textContent,
                   results: !!document.getElementById('viewer-color-results')})""")
    not_read = answers["/api/residue-values"]["not_read"]
    assert answers["/api/residue-values"]["properties"] == []
    assert "analysis/rmsf/options.json was not read" in not_read
    assert said == {"note": f"Per-residue results not read: {not_read}", "results": False}


def test_the_report_page_says_only_that_the_report_was_not_read(browser, tmp_path):
    said, answers = _said(
        browser, _linked(tmp_path), "report",
        "() => document.getElementById('report-empty-reason').textContent",
        """() => [...document.querySelectorAll('#report-empty p')]
                   .filter((p) => !p.hidden).map((p) => p.textContent)""")
    not_read = answers["/api/report"]["not_read"]
    assert not_read.startswith("report/report.md was not read")
    assert said == [not_read]


def test_the_remedy_card_says_its_record_was_not_read(browser, tmp_path):
    root = _stopped(tmp_path)
    _linked_out(root, tmp_path, "manifest.json")
    target = tmp_path / "outside" / "manifest.json"
    target.write_text(target.read_text(encoding="utf-8").replace("SIGTERM", MARK),
                      encoding="utf-8")
    said, answers = _said(
        browser, root, "overview",
        "() => !!document.querySelector('#fixes-card:not([hidden]) .fix')",
        """() => ({text: document.querySelector('#fixes-card .fix').textContent,
                   run: !!document.querySelector('#fixes-card .fix-run')})""")
    assert answers["/api/fixes"]["unavailable"] == ["manifest.json"]
    assert MARK not in said["text"]
    assert said["text"].startswith("The study: manifest.json was not read")
    assert said["run"] is False
