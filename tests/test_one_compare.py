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


def test_a_record_linked_out_of_a_run_compared_is_not_read(variants, tmp_path):
    from tests.test_the_study_pages_read_only_the_study import (
        OUTSIDE, TEXT, _ask, _OPENED, _watch, _WATCHED)

    first, second = _runs(variants)
    outside = tmp_path / "elsewhere"
    _linked_out(first / "analysis" / "rmsd" / "options.json", outside / "options.json",
                OUTSIDE["determined"])
    _watch(outside)
    try:
        status, body = _ask(variants, tmp_path, _compared_route(first, second), public=False)
        opened = sorted(set(_OPENED))
    finally:
        _WATCHED.clear()
    assert opened == []
    assert TEXT not in body.decode("utf-8")
    answer = json.loads(body)
    assert status == 200 and answer["ok"] is True
    rmsd = next(m for m in answer["measures"] if m["analysis"] == "rmsd")
    # Said not read on its side, never as a run that recorded no mean.
    assert rmsd["first"]["not_read"] is True and rmsd["first"]["mean"] is None
    assert "versus" not in rmsd
    assert answer["not_read"].startswith(
        "study, temperature 300 K: analysis/rmsd/options.json was not read")
    assert "rmsd" not in [s["analysis"] for s in answer["series"]]
    # Both recorded it: not drawn because it was not read, never said as a
    # series they do not have in common.
    assert answer["series_not_read"] == ["RMSD"]


def test_a_config_linked_out_is_said_not_read_not_missing(variants, tmp_path):
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
        status, body = _ask(variants, tmp_path, _compared_route(first, second), public=False)
        opened = sorted(set(_OPENED))
    finally:
        _WATCHED.clear()
    assert opened == []
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
    # The two runs' study: its record not read said once, not once a run.
    assert answer["not_read"].count("study: batch_manifest.json was not read") == 1
    # Whether they are repeats of one measurement is then not known: said.
    assert answer["replicas_unknown"].startswith("Whether these two runs are replicas is not known")


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


# ---------------------------------------------------------------------------
# Two replicas: the error between them is the one `batch/aggregate.py` gives.

def test_two_replicas_are_set_against_the_replicas_spread_not_their_own_errors(tmp_path):
    """Three runs differing only by seed, each with its own error ten times
    and more smaller than the spread of their means: the error each run
    printed is too tight (`batch.aggregate`'s calibration), and a difference
    between two of them set against it was called resolved. Set against the
    error `aggregate` gives each replica, it is within the errors."""
    from fastmdxplora.batch.aggregate import aggregate_members
    from fastmdxplora.gui.compare_view import compared
    from tests.test_the_overview_of_a_study_of_runs import _a_study_of_runs

    study = _a_study_of_runs(tmp_path, "simulation.random_seed", [1, 2, 3])
    first, second, _ = _runs(study)
    said = compared(first, second)
    rmsd = next(m for m in said["measures"] if m["analysis"] == "rmsd")
    own = [rmsd["first"]["error"], rmsd["second"]["error"]]
    assert abs(rmsd["versus"]["difference"]) > 10 * (own[0] ** 2 + own[1] ** 2) ** 0.5
    entry = aggregate_members(study)["analyses"]["rmsd"]
    assert entry["ratio_observed_to_predicted"] > 10
    errors = [row["error_as_a_replica"] for row in entry["members"][:2]]
    assert rmsd["versus"]["errors"] == pytest.approx(errors)
    assert rmsd["versus"]["error"] == pytest.approx((errors[0] ** 2 + errors[1] ** 2) ** 0.5)
    assert rmsd["versus"]["resolved"] is False and rmsd["versus"]["error_of"] == "replicas"
    assert said["replicas"]["runs"] == 3


def test_two_runs_of_a_sweep_keep_their_own_errors(variants):
    from fastmdxplora.gui.compare_view import compared

    said = compared(*_runs(variants))
    rmsd = next(m for m in said["measures"] if m["analysis"] == "rmsd")
    assert "error_of" not in rmsd["versus"] and "replicas" not in said


def test_the_replica_error_is_the_spread_or_the_run_s_own_where_larger():
    from fastmdxplora.batch.aggregate import replica_error

    assert replica_error([1.0, 2.0, 3.0], 0.01) == pytest.approx(1.0)
    assert replica_error([1.0, 2.0, 3.0], 5.0) == pytest.approx(5.0)
    assert replica_error([1.0, 2.0, 3.0], None) == pytest.approx(1.0)
    assert replica_error([1.0], 0.5) is None


# ---------------------------------------------------------------------------
# The thermodynamic means and series, as the run's Overview reads them.

def _on_a_clock(run: Path, every_ps: float) -> None:
    """The run's analyses loaded its frames ``every_ps`` apart, as a run
    with a clock records it."""
    path = run / "analysis" / "analysis_manifest.json"
    manifest = json.loads(path.read_text(encoding="utf-8")) if path.exists() else {}
    manifest["load_kwargs"] = {**(manifest.get("load_kwargs") or {}), "saving_interval_ps": every_ps}
    manifest.pop("frame_times_ps", None)
    path.write_text(json.dumps(manifest), encoding="utf-8")


def _energies(run: Path, potential: float) -> None:
    import numpy as np

    noise = np.random.default_rng(1).standard_normal(400)
    lines = ['#"Step","Time (ps)","Potential Energy (kJ/mole)","Temperature (K)"']
    for i in range(400):
        lines.append(f"{(i + 1) * 500},{(i + 1) * 1.0},{potential + 5 * noise[i]:.4f},"
                     f"{300 + 2 * noise[i]:.4f}")
    (run / "simulation" / "energy.csv").write_text("\n".join(lines) + "\n", encoding="utf-8")


def test_the_thermodynamic_means_and_series_are_compared(variants):
    from fastmdxplora.gui.compare_view import compared
    from fastmdxplora.gui.overview_view import overview_payload

    first, second = _runs(variants)
    for run, potential in ((first, -1000.0), (second, -1100.0)):
        _energies(run, potential)
        _on_a_clock(run, 10.0)
    said = compared(first, second)
    rows = {m["analysis"]: m for m in said["measures"]}
    potential = rows["potential_energy"]
    assert (potential["label"], potential["unit"]) == ("Potential energy", "kJ/mol")
    # As the run's Overview says it.
    shown = overview_payload(first)["thermodynamics"]["means"]["potential_energy"]
    assert potential["first"]["mean"] == pytest.approx(shown["mean"])
    assert potential["versus"]["difference"] == pytest.approx(-100.0, abs=2.0)
    assert potential["versus"]["resolved"] is True
    assert list(rows)[-2:] == ["potential_energy", "temperature"]
    drawn = {s["analysis"]: s for s in said["series"]}
    # Production's clock, named as the analyses' series on the card are.
    assert drawn["potential_energy"]["first"]["x_label"] == drawn["rmsd"]["first"]["x_label"] \
        == "Time (ns)"
    assert len(drawn["potential_energy"]["second"]["y"]) == 400


def test_every_series_on_the_card_is_said_on_one_clock_legibly_on_a_phone(variants):
    """RMSD was drawn against "Time (ps)" and, on the same card, the
    energies against "Simulation time (ns)": two clocks for one run. And on
    a phone each chart was scaled down whole, its labels about 5 px."""
    first, second = _runs(variants)
    for run, potential in ((first, -1000.0), (second, -1100.0)):
        _energies(run, potential)
        _on_a_clock(run, 10.0)
    place = urllib.parse.quote(str(first), safe="") + "|" + urllib.parse.quote(str(second), safe="")

    def look(page, url):
        out = {}
        for width in (1440, 390):
            page.set_viewport_size({"width": width, "height": 900})
            page.goto(url + "#compare/" + place, wait_until="domcontentloaded")
            page.wait_for_selector('.compare-series-chart[data-drawn="2"]')
            out[width] = page.evaluate("""() => {
                const charts = [...document.querySelectorAll('.compare-series-block')];
                return {
                    clocks: charts.map((c) => [c.dataset.analysis,
                        c.querySelector('text.runs-overlay-label').textContent]),
                    smallest: Math.min(...[...document.querySelectorAll(
                        '.compare-series-chart text')].map((t) => t.getBoundingClientRect().height)),
                    sideways: document.documentElement.scrollWidth > window.innerWidth,
                };
            }""")
            page.goto(url + "#overview", wait_until="domcontentloaded")
        return out

    said = _browse(variants, look)
    for width in (1440, 390):
        clocks = dict(said[width]["clocks"])
        assert {"rmsd", "potential_energy", "temperature"} <= set(clocks)
        assert set(clocks.values()) == {"Time (ps)"}, clocks
        assert said[width]["sideways"] is False
    # Labels as legible on a phone as at full width: drawn at their size.
    assert said[390]["smallest"] >= 9, said[390]


# ---------------------------------------------------------------------------
# A run whose Manifest was not read is not said to stand anywhere.

def test_a_run_s_state_is_not_known_when_its_manifest_was_not_read(variants, tmp_path):
    from tests.test_the_study_pages_read_only_the_study import _ask, _OPENED, _watch, _WATCHED

    first, second = _runs(variants)
    manifest = first / "manifest.json"
    outside = tmp_path / "elsewhere"
    _linked_out(manifest, outside / "manifest.json",
                json.loads(manifest.read_text(encoding="utf-8")))
    _watch(outside)
    try:
        status, body = _ask(variants, tmp_path, _compared_route(first, second), public=False)
        opened = sorted(set(_OPENED))
    finally:
        _WATCHED.clear()
    assert status == 200 and opened == []
    answer = json.loads(body)
    assert answer["first"]["state"] == "not known: manifest.json was not read"
    assert answer["second"]["state"] == "completed"


# ---------------------------------------------------------------------------
# On a phone, and the page's own figures.

def test_the_page_on_a_phone_breaks_no_word_and_says_the_replicas_error(tmp_path):
    from tests.test_the_overview_of_a_study_of_runs import _a_study_of_runs

    study = _a_study_of_runs(tmp_path, "simulation.random_seed", [1, 2, 3])
    first, second, _ = _runs(study)
    (second / "resolved_config.yml").write_text(
        (second / "resolved_config.yml").read_text(encoding="utf-8")
        + "setup:\n  random_seed: 7\n", encoding="utf-8")
    place = urllib.parse.quote(str(first), safe="") + "|" + urllib.parse.quote(str(second), safe="")

    def look(page, url):
        page.set_viewport_size({"width": 390, "height": 844})
        page.goto(url + "#compare/" + place, wait_until="domcontentloaded")
        page.wait_for_selector('.page[data-page="compare"][data-compared="yes"]')
        return page.evaluate("""() => {
            const words = [...document.querySelectorAll(
                '.compare-means-table th, .compare-settings-table td.compare-setting')];
            return {
                broken: words.filter((n) => n.textContent.split(/[ ._]/).some((w) => {
                    const probe = document.createElement('span');
                    probe.textContent = w;
                    probe.style.whiteSpace = 'nowrap';
                    n.appendChild(probe);
                    const style = getComputedStyle(n);
                    const room = n.clientWidth - parseFloat(style.paddingLeft)
                        - parseFloat(style.paddingRight);
                    const fits = probe.getBoundingClientRect().width <= room + 0.5;
                    probe.remove();
                    return !fits;
                })).map((n) => n.textContent),
                sideways: document.documentElement.scrollWidth > window.innerWidth,
                note: (document.querySelector('.compare-replicas') || {}).textContent || '',
                errorOf: document.querySelector(
                    '.compare-means-table tr[data-analysis="rmsd"] td:last-child').dataset.errorOf,
                small: window.FastMDXCompare.withError
                    ? window.FastMDXCompare.withError(1.5e-6, 2.5e-7, 'nm') : null,
            };
        }""")

    said = _browse(study, look)
    assert said["sideways"] is False
    assert said["broken"] == [], said
    assert said["errorOf"] == "replicas"
    assert said["note"].startswith("These are two of 3 runs that differ only by random seed")
    # An error below a millionth said in its powers of ten, never 0.000000.
    assert said["small"] == "(15.0 \u00b1 2.5)e-7 nm"


def test_ticking_the_second_run_does_not_move_the_table(variants):
    def look(page, url):
        page.goto(url + "#overview", wait_until="domcontentloaded")
        page.wait_for_selector("#overview-runs-table tbody tr")
        top = "() => document.getElementById('overview-runs-table').getBoundingClientRect().top"
        before = page.evaluate(top)
        page.locator(".overview-runs-pick").nth(0).check()
        page.locator(".overview-runs-pick").nth(1).check()
        page.wait_for_selector("#overview-runs-compare:not([hidden])")
        return before, page.evaluate(top)

    before, after = _browse(variants, look)
    assert abs(after - before) < 1, (before, after)


def test_two_runs_of_a_tripeptide_are_played_beside_in_the_viewer(tmp_path):
    """Two seeds of one study of a tripeptide: Play them beside said the
    second "is not a form of this study's protein", since no run of four
    residues alike fits in three."""
    pytest.importorskip("mdtraj")
    import numpy as np

    from tests.test_another_study_is_shown_beside import _peptide
    from tests.test_the_overview_of_a_study_of_runs import _a_study_of_runs

    study = _a_study_of_runs(tmp_path, "simulation.random_seed", [1, 2])
    first, second = _runs(study)
    for seed, run in enumerate((first, second)):
        peptide = _peptide("IVG")
        xyz = np.repeat(peptide.xyz, 10, axis=0)
        xyz += np.random.default_rng(seed).normal(0, 0.01, xyz.shape)
        frames = type(peptide)(xyz.astype(np.float32), peptide.topology)
        frames[0].save_pdb(str(run / "simulation" / "trajectory_topology.pdb"))
        frames.save_dcd(str(run / "simulation" / "production.dcd"))
        (run / "simulation" / "simulation_parameters.json").write_text(json.dumps(
            {"duration_ns_actual": 1.0}), encoding="utf-8")
    place = urllib.parse.quote(str(first), safe="") + "|" + urllib.parse.quote(str(second), safe="")

    def look(page, url):
        page.goto(url + "#compare/" + place, wait_until="domcontentloaded")
        if not page.evaluate("() => !!document.createElement('canvas').getContext('webgl')"):
            pytest.skip("this browser has no WebGL, so the viewer cannot render")
        page.wait_for_selector('.page[data-page="compare"][data-compared="yes"]')
        page.click("#compare-beside")
        page.wait_for_function("() => document.documentElement.dataset.page === 'viewer'")
        page.wait_for_function("() => !document.getElementById('beside-remove').hidden"
                               " || document.getElementById('beside-note').textContent"
                               ".startsWith('Not shown')", timeout=120000)
        return page.evaluate("() => document.getElementById('beside-note').textContent")

    note = _browse(study, look)
    assert not note.startswith("Not shown"), note
    assert "its 3 residues paired" in note, note


@PUBLIC
def test_beyond_loopback_neither_compare_nor_the_timeline_is_offered(variants, public):
    """On a public bind /api/compared and /api/study-timeline are refused,
    yet Compare with... and the Timeline tab were offered, to say "This is
    disabled when the dashboard is bound to a non-loopback address". The
    page is told how it is bound and offers only what will be answered; a
    panel left on the Timeline shows the Log."""
    pytest.importorskip("playwright.sync_api")
    from playwright.sync_api import sync_playwright

    from fastmdxplora.gui.exploration import DashboardRuntime
    from tests.test_a_public_dashboard_answers_only_what_it_lists import _serving

    runtime = DashboardRuntime(workspace_root=variants, exploration_root=variants.parent,
                               active_root=variants)
    with _serving(runtime, allow_control=not public) as url, sync_playwright() as pw:
        browser = pw.chromium.launch()
        page = browser.new_page(viewport={"width": 1440, "height": 900})
        page.set_default_timeout(60000)
        errors: list[str] = []
        page.on("pageerror", lambda error: errors.append(str(error)))
        page.add_init_script("localStorage.setItem('fmx.sideTab', 'timeline')")
        page.goto(url + "#overview", wait_until="domcontentloaded")
        page.wait_for_function("() => window.FastMDXDashboard && "
                               "window.FastMDXDashboard.state.appState.active_run")
        bound = page.evaluate("() => document.body.dataset.beyondLoopback")
        tabs = page.evaluate("""() => [...document.querySelectorAll('.side-tab')]
            .filter((t) => t.offsetParent !== null).map((t) => t.textContent.trim())""")
        pane = page.evaluate("""() => [...document.querySelectorAll('.side-pane')]
            .filter((p) => !p.hidden).map((p) => p.dataset.sidePane)""")
        page.keyboard.press("Control+k")
        page.wait_for_selector("#palette:not([hidden])")
        page.fill("#palette-input", "compare")
        listed = page.evaluate("""() => [...document.querySelectorAll('.palette-item')]
            .map((i) => i.querySelector('.palette-label').textContent)""")
        browser.close()
    assert errors == []
    assert bound == ("true" if public else "false")
    if public:
        assert "Timeline" not in tabs and pane == ["log"]
        assert "Compare with…" not in listed
    else:
        assert "Timeline" in tabs and pane == ["timeline"]
        assert "Compare with…" in listed
