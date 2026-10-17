"""A plotted source frame opens those coordinates, even when it was not sampled.

The browser's bounded playback is deliberately only frames 0, 5 and 10.
Frame 7 remains a source identity through the chart and is loaded as one
additional coordinate frame when selected.
"""

from __future__ import annotations

import json
import struct
import time
import urllib.request
import zlib
from pathlib import Path

import numpy as np
import pytest

md = pytest.importorskip("mdtraj")


def _study(root: Path, *, frame_count: int = 11, residues: int = 1) -> Path:
    simulation = root / "simulation"
    simulation.mkdir(parents=True)
    topology = md.Topology()
    chain = topology.add_chain()
    atoms_per_residue = ("N", "CA", "C", "O", "CB")
    for residue_index in range(residues):
        residue = topology.add_residue("ALA", chain, resSeq=residue_index + 1)
        for name in atoms_per_residue:
            topology.add_atom(name, md.element.carbon, residue)
    n_atoms = topology.n_atoms
    xyz = np.zeros((frame_count, n_atoms, 3), dtype=np.float32)
    xyz[:, :, 0] = np.arange(frame_count, dtype=np.float32)[:, None] * 0.1
    xyz[:, :, 1] = np.arange(n_atoms, dtype=np.float32)[None, :] * 0.025
    xyz[:, :, 2] = 0.1
    trajectory = md.Trajectory(xyz, topology)
    trajectory[0].save_pdb(str(simulation / "trajectory_topology.pdb"))
    trajectory.save_dcd(str(simulation / "production.dcd"))
    (simulation / "live_status.json").write_text(json.dumps(
        {"status": "completed", "stage": "production"}), encoding="utf-8")
    (simulation / "simulation_parameters.json").write_text(json.dumps({
        "parameters": {"integrator": "verlet", "timestep_fs": 2.0},
        "resolved": {"trajectory_interval_steps": 500,
                     "production_steps": frame_count * 500},
    }), encoding="utf-8")
    return root


def _series(root: Path, *, first: int = 0, stride: int = 1, rows: int = 11) -> None:
    analysis = root / "analysis"
    series = analysis / "rmsd"
    series.mkdir(parents=True)
    (analysis / "analysis_manifest.json").write_text(json.dumps({
        "load_kwargs": {"first": first, "stride": stride, "saving_interval_ps": 1.0},
        "resolved": {"trajectory": str(root / "simulation" / "production.dcd")},
        "results": {"rmsd": {"status": "done"}},
    }), encoding="utf-8")
    (series / "rmsd.dat").write_text(
        "# RMSD\n" + "".join(f"{0.01 * frame:.4f}\n" for frame in range(rows)),
        encoding="utf-8")
    (series / "options.json").write_text(json.dumps({"findings": {}}), encoding="utf-8")
    def chunk(kind: bytes, data: bytes) -> bytes:
        return (struct.pack(">I", len(data)) + kind + data
                + struct.pack(">I", zlib.crc32(kind + data) & 0xFFFFFFFF))
    pixels = b"".join(b"\0" + b"\x80\x80\x80" * 65 for _ in range(42))
    (series / "rmsd.png").write_bytes(
        b"\x89PNG\r\n\x1a\n"
        + chunk(b"IHDR", struct.pack(">IIBBBBB", 65, 42, 8, 2, 0, 0, 0))
        + chunk(b"IDAT", zlib.compress(pixels)) + chunk(b"IEND", b""))


def test_exact_frame_routes_validate_source_signature_and_bounds(tmp_path):
    from fastmdxplora.gui.trajectory_frames import exact_frame_info

    root = _study(tmp_path / "study")
    info = exact_frame_info(root, 7, most_frames=3)
    assert info["ok"]
    assert info["frame_indices"] == [0, 5, 10]
    assert info["time_ns"] == pytest.approx(0.008)
    assert "frame=7" in info["exact_url"]
    assert "exact_frame=7" in info["dcd_url"]
    assert exact_frame_info(root, 7, most_frames=3,
                            source_signature="a different trajectory")["ok"] is False
    assert exact_frame_info(root, -1, most_frames=3)["ok"] is False
    assert exact_frame_info(root, 11, most_frames=3)["ok"] is False
    assert exact_frame_info(tmp_path / "missing", 0, most_frames=3)["ok"] is False


def test_exact_xtc_and_dcd_fallback_hold_the_requested_coordinates(tmp_path):
    from fastmdxplora.gui.trajectory_frames import (
        exact_frame_xtc,
        exact_frames_dcd,
        frames_info,
    )

    root = _study(tmp_path / "study")
    info = frames_info(root, most_frames=3)
    exact, details = exact_frame_xtc(root, 7, most_frames=3,
                                     source_signature=info["source_signature"])
    assert details["ok"] and exact
    xtc = tmp_path / "exact.xtc"
    xtc.write_bytes(exact)
    read = md.load_xtc(str(xtc), top=str(root / "simulation" / "frames_topology.pdb"))
    assert read.n_frames == 1
    assert read.xyz[0, 1, 0] == pytest.approx(0.7, abs=0.001)

    combined, fallback = exact_frames_dcd(root, 7, most_frames=3,
                                         source_signature=info["source_signature"])
    assert fallback["ok"] and combined
    dcd = tmp_path / "combined.dcd"
    dcd.write_bytes(combined)
    read = md.load_dcd(str(dcd), top=str(root / "simulation" / "frames_topology.pdb"))
    assert read.n_frames == 4
    assert read.xyz[-1, 1, 0] == pytest.approx(0.7, abs=1e-5)


def test_nonzero_analysis_start_and_stride_keep_source_frame_numbers(tmp_path):
    from fastmdxplora.gui.series import series_payload

    root = _study(tmp_path / "study", frame_count=15)
    _series(root, first=2, stride=2, rows=4)
    data = series_payload(root, "rmsd")
    assert data["frames"] == [4, 6, 8, 10]


@pytest.mark.parametrize("to", ["first", "start"])
def test_exact_frame_uses_the_viewers_backbone_alignment(tmp_path, to):
    from fastmdxplora.gui.trajectory_frames import (
        exact_frame_xtc,
        frames_info,
        superposed_frames,
    )

    root = _study(tmp_path / "study")
    sample = frames_info(root, most_frames=3)
    aligned = superposed_frames(root, "backbone", to=to, smooth=1)
    assert aligned["ok"]
    expected = md.load_dcd(str(root / "simulation" / aligned["file"]),
                           top=str(root / "simulation" / "frames_topology.pdb"))
    body, details = exact_frame_xtc(root, 7, most_frames=3,
                                    source_signature=sample["source_signature"],
                                    on="backbone", to=to)
    assert details["ok"] and body
    output = tmp_path / f"exact-{to}.xtc"
    output.write_bytes(body)
    exact = md.load_xtc(str(output), top=str(root / "simulation" / "frames_topology.pdb"))
    assert exact.n_frames == 1
    assert np.allclose(exact.xyz[0], expected.xyz[0], atol=0.001)


def test_server_exact_frame_response_is_one_frame_and_signature_guarded(tmp_path):
    from fastmdxplora.gui.server import DashboardConfig, start_dashboard_session

    root = _study(tmp_path / "study")
    _series(root)
    session = start_dashboard_session(output=str(root), host="127.0.0.1", port=0,
                                     config=DashboardConfig(max_browser_frames=3))
    try:
        base = session.url.rstrip("/")
        info = json.loads(urllib.request.urlopen(base + "/api/frames-info").read())
        metadata = json.loads(urllib.request.urlopen(
            base + f"/api/exact-frame?frame=7&signature={info['source_signature']}").read())
        assert metadata["ok"] and metadata["frame_indices"] == [0, 5, 10]
        started = time.perf_counter()
        response = urllib.request.urlopen(base + metadata["exact_url"], timeout=30)
        body = response.read()
        elapsed = time.perf_counter() - started
        assert response.status == 200
        assert len(body) == int(response.headers["Content-Length"])
        assert 0 < len(body) < 4096
        assert elapsed >= 0
        stale = urllib.request.urlopen(base + "/api/exact-frame?frame=7&signature=old")
        assert json.loads(stale.read())["ok"] is False
        non_ascii = urllib.request.urlopen(base + "/api/exact-frame?frame=%C2%B2")
        assert json.loads(non_ascii.read())["ok"] is False
    finally:
        session.server.shutdown()


@pytest.mark.parametrize("route, key", [
    ("/api/exact-frame", "frame"),
    ("/structure/exact-frame.xtc", "frame"),
    ("/structure/frames.dcd", "exact_frame"),
])
@pytest.mark.parametrize("invalid", ["oversized", "repeated", "empty"])
def test_exact_frame_routes_refuse_unreadable_or_ambiguous_indices(tmp_path, route, key, invalid):
    from fastmdxplora.gui.server import start_dashboard_session

    session = start_dashboard_session(output=str(_study(tmp_path / "study")),
                                      host="127.0.0.1", port=0)
    try:
        base = session.url.rstrip("/")
        info = json.loads(urllib.request.urlopen(base + "/api/frames-info").read())
        raw = {"oversized": "9" * 5000, "repeated": f"1&{key}=2", "empty": ""}[invalid]
        url = base + route + f"?{key}={raw}&signature={info['source_signature']}"
        if route.startswith("/api/"):
            with urllib.request.urlopen(url) as response:
                assert json.loads(response.read())["ok"] is False
        else:
            with pytest.raises(urllib.error.HTTPError) as caught:
                urllib.request.urlopen(url)
            assert caught.value.code == 400
        with urllib.request.urlopen(base + "/api/exact-frame?frame=1") as response:
            assert json.loads(response.read())["ok"] is True
    finally:
        session.server.shutdown()


def test_a_chart_click_shows_exact_frame_coordinates_and_playback_order(tmp_path, monkeypatch):
    pytest.importorskip("playwright.sync_api")
    from playwright.sync_api import sync_playwright

    from fastmdxplora.gui import trajectory_frames
    from fastmdxplora.gui.server import DashboardConfig, start_dashboard_session

    root = _study(tmp_path / "study")
    _series(root)
    session = start_dashboard_session(output=str(root), host="127.0.0.1", port=0,
                                     config=DashboardConfig(max_browser_frames=3))
    try:
        with sync_playwright() as pw:
            browser = pw.chromium.launch(args=["--enable-unsafe-swiftshader"])
            page = browser.new_page(viewport={"width": 1440, "height": 900})
            page.add_init_script(
                "try { localStorage.setItem('fmx.viewerPlaybackOpen', '1'); } catch (e) {}")
            errors: list[str] = []
            page.on("pageerror", lambda error: errors.append(str(error)))
            page.goto(session.url + "#analysis", wait_until="domcontentloaded")
            if not page.evaluate("() => !!document.createElement('canvas').getContext('webgl')"):
                pytest.skip("this browser has no WebGL, so the viewer cannot render")
            page.wait_for_selector('.series-chart[data-analysis="rmsd"] svg')
            payload = json.loads(urllib.request.urlopen(
                session.url.rstrip("/") + "/api/frames-info").read())
            assert payload["frame_indices"] == [0, 5, 10]
            chart = page.locator('.series-chart[data-analysis="rmsd"] svg')
            chart.scroll_into_view_if_needed()
            box = chart.bounding_box()
            x = box["x"] + 64 + (box["width"] - 64 - 28) * 0.7
            page.mouse.click(x, box["y"] + box["height"] / 2)
            page.wait_for_function(
                "() => window.FastMDXMoleculeViewer?.STATE.exactFrame?.frame === 7",
                timeout=60000)
            page.wait_for_function(
                "() => window.FastMDXMoleculeViewer.STATE.engine"
                ".find({atom: 'CA'})[0]?.x > 6.98", timeout=90000)
            got = page.evaluate(
                "() => ({frame: document.getElementById('traj-current').textContent,"
                " slider: document.getElementById('traj-slider').value,"
                " sources: window.FastMDXMoleculeViewer.STATE.playbackFrameIndices,"
                " engine: window.FastMDXMoleculeViewer.STATE.playbackEngineIndices,"
                " x: window.FastMDXMoleculeViewer.STATE.engine.find({atom:'CA'})[0].x})"
            )
            assert got["frame"] == "3"
            assert got["sources"] == [0, 5, 7, 10]
            assert got["engine"] == [0, 1, 3, 2]
            assert got["x"] == pytest.approx(7.0, abs=0.02)
            # The next slider step is source frame 10 and the previous one is 5.
            page.locator("#traj-slider").evaluate(
                "el => { el.value = '3'; el.dispatchEvent(new Event('input')); }"
            )
            page.wait_for_function(
                "() => document.getElementById('traj-current').textContent === '4'"
            )
            page.locator("#traj-slider").evaluate(
                "el => { el.value = '1'; el.dispatchEvent(new Event('input')); }"
            )
            page.wait_for_function(
                "() => document.getElementById('traj-current').textContent === '2'"
            )

            def select_chart_frame(frame: int) -> None:
                page.evaluate("() => window.FastMDXDashboard.navigate('analysis')")
                page.wait_for_function("() => document.documentElement.dataset.page === 'analysis'")
                page.wait_for_selector('.series-chart[data-analysis="rmsd"] svg', state="visible")
                svg = page.locator('.series-chart[data-analysis="rmsd"] svg')
                svg.scroll_into_view_if_needed()
                rect = svg.bounding_box()
                point_x = rect["x"] + 64 + (rect["width"] - 64 - 28) * frame / 10
                page.mouse.click(point_x, rect["y"] + rect["height"] / 2)

            select_chart_frame(0)
            page.wait_for_function(
                "() => document.getElementById('traj-current').textContent === '1'"
            )
            page.wait_for_function("() => window.FastMDXMoleculeViewer.STATE.engine"
                ".find({atom:'CA'})[0]?.x < 0.02")
            select_chart_frame(10)
            page.wait_for_function(
                "() => document.getElementById('traj-current').textContent === '4'"
            )
            page.wait_for_function("() => window.FastMDXMoleculeViewer.STATE.engine"
                ".find({atom:'CA'})[0]?.x > 9.98")

            original_info = trajectory_frames.exact_frame_info
            delayed_frames: set[int] = set()
            def exact_info_with_delay(*args, **kwargs):
                frame = args[1] if len(args) > 1 else kwargs.get("frame")
                if frame in delayed_frames:
                    time.sleep(0.7)
                return original_info(*args, **kwargs)
            monkeypatch.setattr(trajectory_frames, "exact_frame_info", exact_info_with_delay)
            delayed_frames.add(8)
            page.evaluate("() => window.FastMDXDashboard.navigate('analysis')")
            page.wait_for_function("() => document.documentElement.dataset.page === 'analysis'")
            page.wait_for_selector('.series-chart[data-analysis="rmsd"] svg', state="visible")
            svg = page.locator('.series-chart[data-analysis="rmsd"] svg')
            svg.scroll_into_view_if_needed()
            rect = svg.bounding_box()
            with page.expect_request("**/api/exact-frame?frame=8*"):
                page.mouse.click(rect["x"] + 64 + (rect["width"] - 64 - 28) * 0.8,
                                 rect["y"] + rect["height"] / 2)
            page.evaluate("() => window.FastMDXDashboard.navigate('analysis')")
            page.wait_for_function("() => document.documentElement.dataset.page === 'analysis'")
            page.wait_for_selector('.series-chart[data-analysis="rmsd"] svg', state="visible")
            svg = page.locator('.series-chart[data-analysis="rmsd"] svg')
            rect = svg.bounding_box()
            page.mouse.click(rect["x"] + 64 + (rect["width"] - 64 - 28) * 0.9,
                             rect["y"] + rect["height"] / 2)
            page.wait_for_function(
                "() => window.FastMDXMoleculeViewer.STATE.exactFrame?.frame === 9", timeout=60000
            )
            page.wait_for_timeout(1000)
            assert page.evaluate("() => window.FastMDXMoleculeViewer.STATE.exactFrame.frame") == 9
            assert page.evaluate(
                "() => window.FastMDXMoleculeViewer.STATE.engine.find({atom:'CA'})[0].x"
            ) == pytest.approx(9.0, abs=0.02)

            original_xtc = trajectory_frames.exact_frame_xtc
            def interrupted_frame(output_dir, frame, **kwargs):
                if frame == 6:
                    return None, {"ok": False, "reason": "simulated interrupted request"}
                return original_xtc(output_dir, frame, **kwargs)
            monkeypatch.setattr(trajectory_frames, "exact_frame_xtc", interrupted_frame)
            select_chart_frame(6)
            page.wait_for_function("() => document.getElementById('sr-live').textContent"
                ".includes('Frame 6 could not be shown')", timeout=30000)
            assert page.evaluate("() => window.FastMDXMoleculeViewer.STATE.exactFrame.frame") == 9

            delayed_frames.add(6)
            monkeypatch.setattr(trajectory_frames, "exact_frame_xtc", original_xtc)
            page.evaluate("() => window.FastMDXDashboard.navigate('analysis')")
            page.wait_for_function("() => document.documentElement.dataset.page === 'analysis'")
            page.wait_for_selector('.series-chart[data-analysis="rmsd"] svg', state="visible")
            svg = page.locator('.series-chart[data-analysis="rmsd"] svg')
            svg.scroll_into_view_if_needed()
            rect = svg.bounding_box()
            with page.expect_request("**/api/exact-frame?frame=6*"):
                page.mouse.click(rect["x"] + 64 + (rect["width"] - 64 - 28) * 0.6,
                                 rect["y"] + rect["height"] / 2)
            page.evaluate("() => window.dispatchEvent(new CustomEvent('dashboard:run-changed',"
                "{detail: {first: false}}))")
            page.wait_for_timeout(1100)
            assert page.evaluate("() => window.FastMDXMoleculeViewer.STATE.exactFrame") is None
            assert errors == []
            browser.close()
    finally:
        session.server.shutdown()


def test_figure_chart_supports_keyboard_navigation_and_ignores_clicks_outside_its_plot(tmp_path):
    pytest.importorskip("playwright.sync_api")
    from playwright.sync_api import sync_playwright

    from fastmdxplora.gui.server import DashboardConfig, start_dashboard_session
    from tests.test_an_analysis_is_drawn_from_its_numbers import _figure

    root = _study(tmp_path / "study")
    _series(root)
    folder = root / "analysis" / "cluster"
    folder.mkdir(parents=True)
    (folder / "cluster_kmeans.dat").write_text(
        "frame,cluster\n" + "".join(f"{frame},{frame % 2}\n" for frame in range(11)),
        encoding="utf-8")
    _figure(folder / "cluster_kmeans.png")
    manifest_path = root / "analysis" / "analysis_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["results"]["cluster_kmeans"] = {"status": "done"}
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")

    session = start_dashboard_session(output=str(root), host="127.0.0.1", port=0,
                                     config=DashboardConfig(max_browser_frames=3))
    try:
        with sync_playwright() as pw:
            browser = pw.chromium.launch(args=["--enable-unsafe-swiftshader"])
            page = browser.new_page(viewport={"width": 1440, "height": 900})
            page.add_init_script(
                "try { localStorage.setItem('fmx.viewerPlaybackOpen', '1'); } catch (e) {}")
            errors: list[str] = []
            page.on("pageerror", lambda error: errors.append(str(error)))
            page.goto(session.url + "#analysis", wait_until="domcontentloaded")
            if not page.evaluate("() => !!document.createElement('canvas').getContext('webgl')"):
                pytest.skip("this browser has no WebGL, so the viewer cannot render")
            chart = page.locator('.figure-chart[data-figure-kind="labels"] svg')
            chart.wait_for()
            chart.scroll_into_view_if_needed()
            box = chart.bounding_box()
            point_x = box["x"] + 62 + (box["width"] - 82) * 0.7
            point_y = box["y"] + box["height"] * 0.5
            page.mouse.click(point_x, point_y)
            page.wait_for_function(
                "() => window.FastMDXMoleculeViewer?.STATE.exactFrame?.frame === 7",
                timeout=60000)
            page.wait_for_function(
                "() => window.FastMDXMoleculeViewer.STATE.engine.find({atom:'CA'})[0]?.x > 6.98",
                timeout=90000)

            # Reopen the analysis chart after the exact-frame action so its
            # keyboard state is initialized against the updated playback map.
            page.evaluate("() => window.FastMDXDashboard.navigate('viewer')")
            page.wait_for_function("() => document.documentElement.dataset.page === 'viewer'")
            page.evaluate("() => window.FastMDXDashboard.navigate('analysis')")
            page.wait_for_function("() => document.documentElement.dataset.page === 'analysis'")
            chart = page.locator('.figure-chart[data-figure-kind="labels"] svg')
            chart.wait_for(state="visible")
            page.wait_for_function("() => document.querySelector("
                                   "'.figure-chart[data-figure-kind=labels] svg')"
                                   "?.getAttribute('role') === 'application'")
            chart.scroll_into_view_if_needed()
            chart.focus()
            chart.press("End")
            chart.press("Enter")
            page.wait_for_function(
                "() => { const s = window.FastMDXMoleculeViewer.STATE;"
                " const i = Number(document.getElementById('traj-slider').value);"
                " return s.playbackFrameIndices[i] === 10; }",
                timeout=60000)
            page.wait_for_function(
                "() => window.FastMDXMoleculeViewer.STATE.engine.find({atom:'CA'})[0]?.x > 9.98",
                timeout=90000)

            # Frame navigation can put the viewer page in front. Return to
            # Analysis and wait for the replacement chart to be visible before
            # checking that a click outside its plot leaves the frame alone.
            page.evaluate("() => window.FastMDXDashboard.navigate('analysis')")
            page.wait_for_function("() => document.documentElement.dataset.page === 'analysis'")
            chart = page.locator('.figure-chart[data-figure-kind="labels"] svg')
            chart.wait_for(state="visible")
            chart.scroll_into_view_if_needed()
            box = chart.bounding_box()
            page.mouse.move(box["x"] + 62 + (box["width"] - 82) * 0.2, point_y)
            page.mouse.click(box["x"] + 30, point_y)
            page.wait_for_timeout(700)
            stayed = page.evaluate(
                "() => { const s = window.FastMDXMoleculeViewer.STATE;"
                " const i = Number(document.getElementById('traj-slider').value);"
                " return {frame: s.playbackFrameIndices[i],"
                " x: s.engine.find({atom:'CA'})[0].x}; }"
            )
            browser.close()
    finally:
        session.server.shutdown()

    assert errors == []
    assert stayed["frame"] == 10
    assert stayed["x"] == pytest.approx(10.0, abs=0.02)


def test_a_larger_exact_frame_remains_bounded(tmp_path):
    from fastmdxplora.gui.trajectory_frames import exact_frame_xtc, frames_info

    root = _study(tmp_path / "study", frame_count=101, residues=2000)
    info = frames_info(root, most_frames=3)
    started = time.perf_counter()
    body, details = exact_frame_xtc(root, 73, most_frames=3,
                                    source_signature=info["source_signature"])
    elapsed = time.perf_counter() - started
    assert details["ok"] and body
    assert info["frame_indices"] == [0, 50, 100]
    assert len(body) < (info["n_atoms"] * 12)
    output = tmp_path / "exact.xtc"
    output.write_bytes(body)
    exact = md.load_xtc(str(output), top=str(root / "simulation" / "frames_topology.pdb"))
    assert exact.n_frames == 1
    assert exact.xyz[0, 1, 0] * 10 == pytest.approx(73.0, abs=0.01)
    assert elapsed >= 0
