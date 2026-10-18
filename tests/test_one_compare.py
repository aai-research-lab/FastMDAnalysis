"""One Compare, for two studies or two runs, opened from wherever two are.

Comparison was in three places with three UIs: two studies chosen in All
studies drew their settings and means under the cards, a study of several
runs drew its own table and overlay on its Analysis page, and the Viewer
played another study beside the one open, chosen there from a list without
runs. Compare is now one page: the two named (a run by its study and by
what sets it apart), the settings that differ, the means with their
difference marked only where it is resolved, each series both recorded
overlaid, and the Viewer playing the second beside the first. It is opened
from All studies, from two runs chosen on a study's Overview and from the
palette (Compare with...), and its address names the two, so Back and a
reload keep them.
"""

from __future__ import annotations

import json
import shutil
import urllib.parse
import urllib.request
from pathlib import Path

import pytest

from tests.test_the_overview_of_a_study_of_runs import _a_study_of_runs


@pytest.fixture()
def variants(tmp_path) -> Path:
    return _a_study_of_runs(tmp_path, "simulation.temperature_K", [300, 310])


def _runs(study: Path) -> list[Path]:
    return sorted((study / "runs").iterdir())


def test_two_runs_are_named_by_their_study_and_compared_with_their_series(variants):
    from fastmdxplora.gui.compare_view import compared

    first, second = _runs(variants)
    said = compared(first, second)
    assert said["ok"] is True
    assert said["first"]["said"] == "study, temperature 300 K"
    assert said["second"]["said"] == "study, temperature 310 K"
    assert said["first"]["run_of"]["path"] == str(variants)
    # The runs here are copies of one finished study: their records differ
    # only in the means each recorded.
    assert said["settings"] == []
    rmsd = next(m for m in said["measures"] if m["analysis"] == "rmsd")
    assert rmsd["versus"]["difference"] == pytest.approx(0.01)
    assert [s["analysis"] for s in said["series"]] == ["rmsd", "rg"]
    for s in said["series"]:
        assert s["first"]["x_label"] == s["second"]["x_label"]
        assert len(s["first"]["x"]) == len(s["first"]["y"]) > 1


def test_a_study_is_named_by_its_folder_and_a_pair_that_cannot_be_is_said(variants, tmp_path):
    from fastmdxplora.gui.compare_view import compared

    other = tmp_path / "other"
    shutil.copytree(_runs(variants)[0], other)
    said = compared(other, _runs(variants)[1])
    assert said["first"]["said"] == "other" and said["first"]["run_of"] is None
    refused = compared(tmp_path / "nowhere", other)
    assert refused["ok"] is False and refused["reason"]


def test_the_route_is_answered_here_and_not_beyond_loopback(variants):
    from fastmdxplora.gui.server import GETS_ANSWERED_BEYOND_LOOPBACK, start_dashboard_session

    assert "/api/compared" not in GETS_ANSWERED_BEYOND_LOOPBACK
    first, second = _runs(variants)
    session = start_dashboard_session(output=str(variants), host="127.0.0.1", port=0)
    try:
        answer = json.loads(urllib.request.urlopen(
            session.url.rstrip("/") + "/api/compared?" + urllib.parse.urlencode(
                {"a": str(first), "b": str(second)}), timeout=30).read())
    finally:
        session.server.shutdown()
    assert answer["ok"] is True and len(answer["series"]) == 2


PAGE = """() => {
    const page = document.querySelector('.page[data-page="compare"]');
    return {
        page: document.documentElement.dataset.page, compared: page.dataset.compared || '',
        said: document.getElementById('compare-said').textContent,
        heads: [...document.querySelectorAll('.compare-head-name')].map((n) => n.textContent),
        settings: document.querySelectorAll('.compare-settings-table tr').length,
        means: [...document.querySelectorAll('.compare-means-table tr[data-analysis]')]
            .map((r) => r.dataset.analysis),
        lines: [...document.querySelectorAll('.compare-series-chart')]
            .map((c) => c.querySelectorAll('polyline').length),
    };
}"""


def _browse(study: Path, look):
    pytest.importorskip("playwright.sync_api")
    from playwright.sync_api import sync_playwright

    from fastmdxplora.gui.server import start_dashboard_session

    session = start_dashboard_session(output=str(study), host="127.0.0.1", port=0)
    try:
        with sync_playwright() as pw:
            browser = pw.chromium.launch(args=["--enable-unsafe-swiftshader"])
            page = browser.new_page(viewport={"width": 1440, "height": 900})
            page.set_default_timeout(60000)
            errors: list[str] = []
            page.on("pageerror", lambda error: errors.append(str(error)))
            said = look(page, session.url)
            browser.close()
    finally:
        session.server.shutdown()
    assert errors == []
    return said


def test_two_runs_chosen_on_the_overview_open_compare_and_back_returns(variants):
    def look(page, url):
        page.goto(url + "#overview", wait_until="domcontentloaded")
        page.wait_for_selector("#overview-runs-table tbody tr")
        page.locator(".overview-runs-pick").nth(0).check()
        offered_one = page.is_visible("#overview-runs-compare")
        page.locator(".overview-runs-pick").nth(1).check()
        page.click("#overview-runs-compare")
        page.wait_for_selector('.page[data-page="compare"][data-compared="yes"]')
        page.wait_for_function(
            "() => document.querySelectorAll('.compare-series-chart polyline').length === 4")
        compared = page.evaluate(PAGE)
        address = page.evaluate("() => location.hash")
        page.reload(wait_until="domcontentloaded")
        page.wait_for_selector('.page[data-page="compare"][data-compared="yes"]')
        reloaded = page.evaluate(PAGE)
        page.click("#compare-swap")
        page.wait_for_function("() => document.getElementById('compare-said').textContent"
                               ".startsWith('study, temperature 310 K')")
        page.go_back()
        page.wait_for_function("() => document.getElementById('compare-said').textContent"
                               ".startsWith('study, temperature 300 K')")
        page.go_back()
        page.wait_for_function("() => document.documentElement.dataset.page === 'overview'")
        return offered_one, compared, address, reloaded

    offered_one, compared, address, reloaded = _browse(variants, look)
    assert offered_one is False
    assert compared["page"] == "compare" and compared["compared"] == "yes"
    assert compared["heads"] == ["study, temperature 300 K", "study, temperature 310 K"]
    assert compared["said"] == "study, temperature 300 K and study, temperature 310 K"
    assert compared["means"][:2] == ["rmsd", "rg"]
    assert compared["lines"] == [2, 2]
    assert address.startswith("#compare/") and "|" in address
    assert reloaded["heads"] == compared["heads"], reloaded


def test_compare_with_in_the_palette_lists_the_other_runs_and_studies(variants):
    def look(page, url):
        page.goto(url + "#overview", wait_until="domcontentloaded")
        page.wait_for_function("() => window.FastMDXDashboard && "
                               "window.FastMDXDashboard.state.appState.active_run")
        page.keyboard.press("Control+k")
        page.wait_for_selector("#palette:not([hidden])")
        page.fill("#palette-input", "compare with")
        page.keyboard.press("Enter")
        page.wait_for_function("() => document.getElementById('palette-input').placeholder"
                               ".startsWith('A study, or a run, to compare')")
        page.wait_for_selector(".palette-item")
        listed = page.evaluate("""() => [...document.querySelectorAll('.palette-item')]
            .map((i) => i.querySelector('.palette-label').textContent)""")
        page.fill("#palette-input", "temperature 300")
        page.keyboard.press("Enter")
        page.wait_for_selector('.page[data-page="compare"][data-compared="yes"]')
        return listed, page.evaluate(PAGE)

    listed, compared = _browse(variants, look)
    assert "temperature 300 K" in listed and "temperature 310 K" in listed
    assert compared["heads"] == ["study", "study, temperature 300 K"], compared


def test_the_viewer_plays_the_second_beside_the_first(tmp_path):
    from tests.test_the_cartoon_is_dssp_of_each_frame import _helical_study

    first = _helical_study(tmp_path / "work" / "alpha", analyse=False)
    (first / "resolved_config.yml").write_text(
        "systems:\n  - system: 1HHO\nsimulation:\n  temperature_K: 300\n", encoding="utf-8")
    second = tmp_path / "work" / "beta"
    shutil.copytree(first, second)
    (second / "resolved_config.yml").write_text(
        "systems:\n  - system: 1HHO\nsimulation:\n  temperature_K: 310\n", encoding="utf-8")
    place = urllib.parse.quote(str(first), safe="") + "|" + urllib.parse.quote(str(second), safe="")

    def look(page, url):
        page.goto(url + "#compare/" + place, wait_until="domcontentloaded")
        if not page.evaluate("() => !!document.createElement('canvas').getContext('webgl')"):
            pytest.skip("this browser has no WebGL, so the viewer cannot render")
        page.wait_for_selector('.page[data-page="compare"][data-compared="yes"]')
        page.click("#compare-beside")
        page.wait_for_function("() => document.documentElement.dataset.page === 'viewer'")
        page.wait_for_function("() => !document.getElementById('beside-remove').hidden",
                               timeout=120000)
        return page.evaluate("""() => [document.getElementById('beside-study').value,
            document.getElementById('beside-note').textContent,
            window.FastMDXViewerRail.chosen]""")

    value, note, chosen = _browse(first, look)
    assert value == str(second)
    assert note.startswith("beta,") and "paired" in note
    assert chosen == "side-beside"


# ---------------------------------------------------------------------------
# Each of the two read held to itself: a record linked out of it is not read.

PUBLIC = pytest.mark.parametrize("public", [False, True], ids=["loopback", "public"])


def _compared_route(first: Path, second: Path) -> str:
    return "/api/compared?" + urllib.parse.urlencode({"a": str(first), "b": str(second)})


def _linked_out(record: Path, outside: Path, content: dict) -> None:
    from tests.test_the_study_pages_read_only_the_study import _link

    outside.parent.mkdir(parents=True, exist_ok=True)
    outside.write_text(json.dumps(content), encoding="utf-8")
    record.unlink()
    _link(record, outside)


@PUBLIC
def test_a_record_linked_out_of_a_run_compared_is_not_read(variants, tmp_path, public):
    from tests.test_the_study_pages_read_only_the_study import (
        OUTSIDE, TEXT, _ask, _OPENED, _watch, _WATCHED)

    first, second = _runs(variants)
    outside = tmp_path / "elsewhere"
    _linked_out(first / "analysis" / "rmsd" / "options.json", outside / "options.json",
                OUTSIDE["determined"])
    _watch(outside)
    try:
        status, body = _ask(variants, tmp_path, _compared_route(first, second), public=public)
        opened = sorted(set(_OPENED))
    finally:
        _WATCHED.clear()
    assert opened == []
    assert TEXT not in body.decode("utf-8")
    if public:
        # Not answered beyond loopback at all.
        assert status == 403
        return
    answer = json.loads(body)
    assert status == 200 and answer["ok"] is True
    rmsd = next(m for m in answer["measures"] if m["analysis"] == "rmsd")
    # Said not read on its side, never as a run that recorded no mean.
    assert rmsd["first"]["not_read"] is True and rmsd["first"]["mean"] is None
    assert "versus" not in rmsd
    assert answer["not_read"].startswith(
        "study, temperature 300 K: analysis/rmsd/options.json was not read")
    assert "rmsd" not in [s["analysis"] for s in answer["series"]]


@PUBLIC
def test_a_config_linked_out_is_said_not_read_not_missing(variants, tmp_path, public):
    from tests.test_the_study_pages_read_only_the_study import _ask, _OPENED, _watch, _WATCHED

    first, second = _runs(variants)
    outside = tmp_path / "elsewhere"
    config = first / "resolved_config.yml"
    (outside).mkdir()
    (outside / "resolved_config.yml").write_text(config.read_text(encoding="utf-8"),
                                                 encoding="utf-8")
    config.unlink()
    from tests.test_the_study_pages_read_only_the_study import _link

    _link(config, outside / "resolved_config.yml")
    _watch(outside)
    try:
        status, body = _ask(variants, tmp_path, _compared_route(first, second), public=public)
        opened = sorted(set(_OPENED))
    finally:
        _WATCHED.clear()
    assert opened == []
    if public:
        assert status == 403
        return
    answer = json.loads(body)
    assert answer["ok"] is False
    assert answer["reason"].startswith(f"{first.name}: resolved_config.yml was not read")
    assert "holds no resolved_config.yml" not in answer["reason"]


def test_a_run_s_study_record_linked_out_does_not_name_it(variants, tmp_path):
    """Where a run is is read from its study's records, held to that study."""
    from tests.test_the_study_pages_read_only_the_study import _ask, _OPENED, _watch, _WATCHED

    first, second = _runs(variants)
    manifest = variants / "batch_manifest.json"
    outside = tmp_path / "elsewhere"
    _linked_out(manifest, outside / "batch_manifest.json",
                json.loads(manifest.read_text(encoding="utf-8")))
    _watch(outside)
    try:
        status, body = _ask(variants, tmp_path, _compared_route(first, second), public=False)
        opened = sorted(set(_OPENED))
    finally:
        _WATCHED.clear()
    assert status == 200 and opened == []
    answer = json.loads(body)
    assert answer["first"]["said"] == first.name and answer["first"]["run_of"] is None
    assert "study: batch_manifest.json was not read" in answer["not_read"]


def test_the_page_says_a_record_not_read(variants, tmp_path):
    from tests.test_the_study_pages_read_only_the_study import OUTSIDE

    first, second = _runs(variants)
    _linked_out(first / "analysis" / "rmsd" / "options.json",
                tmp_path / "elsewhere" / "options.json", OUTSIDE["determined"])
    place = urllib.parse.quote(str(first), safe="") + "|" + urllib.parse.quote(str(second), safe="")

    def look(page, url):
        page.goto(url + "#compare/" + place, wait_until="domcontentloaded")
        page.wait_for_selector('.page[data-page="compare"][data-compared="yes"]')
        return page.evaluate("""() => ({
            line: document.getElementById('compare-not-read').textContent,
            shown: !document.getElementById('compare-not-read').hidden,
            rmsd: [...document.querySelectorAll('.compare-means-table tr[data-analysis="rmsd"] td')]
                .map((td) => td.textContent)})""")

    said = _browse(variants, look)
    assert said["shown"] and "analysis/rmsd/options.json was not read" in said["line"]
    assert said["rmsd"][1] == "not read" and said["rmsd"][3] == "not compared"
