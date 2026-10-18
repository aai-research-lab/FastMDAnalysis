"""The Viewer's playback as a timeline, always under the molecule.

It was folded under the canvas, a slider and a time, with one series
beside it when opened. It is now always there while there are frames to
play: the transport, a ruler of the frames on production's clock with what
happened on them marked (each cluster's medoid, each judgement of a study run
until it is determined), the averaging's start shaded, and the study's series
as tracks moving with the frame; a click on the ruler, a mark or a track
shows that frame. The tool beside the molecule is a panel of its own, as
tall as the rail, not a fold forced open.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest


def _medoids(root: Path, rows: list[tuple[int, int, float]]) -> None:
    folder = root / "analysis" / "cluster"
    folder.mkdir(parents=True, exist_ok=True)
    with (folder / "cluster_kmeans_populations.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["cluster", "frames", "fraction", "medoid_frame", "medoid_time_ns"])
        for cluster, frame, share in rows:
            writer.writerow([cluster, 3, share, frame, 0.002 * frame])


def test_the_marks_are_read_from_the_records(tmp_path):
    from fastmdxplora.gui.timeline_view import timeline_payload

    _medoids(tmp_path, [(0, 2, 0.5), (1, 4, 0.5)])
    (tmp_path / "stopping.json").write_text(json.dumps({
        "targets": [{"analysis": "rmsd", "standard_error": 0.01}],
        "rounds": [{"production_ns": 0.004, "decision": "extend"},
                   {"production_ns": 0.008, "decision": "determined"}],
        "outcome": "determined"}), encoding="utf-8")
    marks = timeline_payload(tmp_path)["marks"]
    medoids = [m for m in marks if m["kind"] == "medoid"]
    checks = [m for m in marks if m["kind"] == "check"]
    assert [m["frame"] for m in medoids] == [2, 4]
    assert medoids[0]["label"] == "The medoid of kmeans cluster 0, 50% of the frames"
    assert [m["time_ns"] for m in checks] == [0.004, 0.008]
    assert checks[1]["label"] == "Judged after round 2: determined"
    assert timeline_payload(tmp_path / "nothing") == {"ok": True, "marks": []}


@pytest.fixture(scope="module")
def study(tmp_path_factory) -> Path:
    from tests.test_the_cartoon_is_dssp_of_each_frame import _helical_study

    root = _helical_study(tmp_path_factory.mktemp("timeline") / "study", analyse=False)
    from fastmdxplora.analysis import AnalysisOrchestrator

    AnalysisOrchestrator(str(root / "simulation" / "production.dcd"),
                         str(root / "simulation" / "trajectory_topology.pdb"),
                         output_dir=str(root / "analysis")).run(include=["rmsd", "rg"])
    _medoids(root, [(0, 3, 0.6), (1, 5, 0.4)])
    return root


def test_the_timeline_is_under_the_molecule_and_shows_a_frame_when_asked(study):
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
            page.goto(session.url + "#viewer", wait_until="domcontentloaded")
            if not page.evaluate("() => !!document.createElement('canvas').getContext('webgl')"):
                pytest.skip("this browser has no WebGL, so the viewer cannot render")
            # Shown without opening anything, and before any frame is shown.
            page.wait_for_selector("#timeline-ruler:not([hidden]) svg")
            page.wait_for_function(
                "() => document.querySelectorAll('#frame-series .frame-series-line').length === 2")
            first = page.evaluate("""() => ({
                fold: !!document.querySelector('details#viewer-under-fold, #viewer-timeline > summary'),
                tracks: [...document.querySelectorAll('.frame-track')].map((t) => t.dataset.analysis),
                marks: [...document.querySelectorAll('.timeline-mark')].map((m) =>
                    m.getAttribute('aria-label')),
            })""")
            # A medoid's mark shows its frame.
            page.locator('.timeline-mark[data-kind="medoid"]').nth(1).dispatch_event("pointerdown")
            page.wait_for_function("() => document.getElementById('timeline-ruler')"
                                   ".dataset.frame === '5'")
            on_mark = page.evaluate("() => [document.getElementById('traj-slider').value,"
                                    " [...document.querySelectorAll('.frame-track-value')]"
                                    ".map((v) => v.textContent)]")
            # The ruler's left end shows the first frame; a track's, too.
            page.locator("#timeline-ruler svg").click(position={"x": 49, "y": 14})
            page.wait_for_function("() => document.getElementById('traj-slider').value === '0'")
            track = page.locator(".frame-track").nth(1).locator("svg")
            width = track.evaluate("(svg) => svg.getBoundingClientRect().width")
            track.click(position={"x": width - 11, "y": 25})
            page.wait_for_function("() => Number(document.getElementById('traj-slider').value) >= 4")
            # A track is taken away and added back.
            page.click('.frame-track[data-analysis="rg"] .frame-track-drop')
            page.wait_for_function("() => document.querySelectorAll('.frame-track').length === 1")
            page.select_option("#frame-series-pick", "rg")
            page.wait_for_function(
                "() => document.querySelectorAll('#frame-series .frame-series-line').length === 2")
            # The tool beside the molecule is as tall as the rail.
            panel = page.evaluate("""() => {
                const chosen = document.querySelector('.viewer-side .side-section.is-chosen');
                const rail = document.getElementById('viewer-rail');
                const head = chosen.querySelector(':scope > summary');
                return [Math.round(chosen.getBoundingClientRect().height),
                        Math.round(rail.getBoundingClientRect().height), head.tabIndex];
            }""")
            browser.close()
    finally:
        session.server.shutdown()
    assert errors == []
    assert first["fold"] is False
    assert first["tracks"] == ["rmsd", "rg"]
    assert first["marks"] == ["The medoid of kmeans cluster 0, 60% of the frames",
                              "The medoid of kmeans cluster 1, 40% of the frames"]
    assert on_mark[0] == "5" and all(v.endswith(" nm") for v in on_mark[1])
    assert panel[0] >= panel[1] - 2 and panel[2] == -1
