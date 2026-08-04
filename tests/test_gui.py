"""Tests for the local live dashboard.

These tests assert:

  - endpoint contracts (URLs, response keys) the frontend binds to,
  - failure isolation (traversal, missing files, empty run),
  - the redesigned HTML/CSS/JS assets exist on disk and contain the
    expected AAI-branded module structure.
"""

from __future__ import annotations

import json
import socket
from pathlib import Path
from urllib.error import HTTPError
from urllib.request import urlopen
from types import SimpleNamespace

import pytest

from fastmdxplora.cli.main import _build_parser, _enable_dashboard_telemetry
from fastmdxplora.gui import protein_preview
from fastmdxplora.gui.ligand_detection import detect_ligands, normalise_ligand_resname
from fastmdxplora.gui.live_frames import (
    dashboard_display_pdb,
    read_live_frame_history,
    write_live_frame,
    write_openmm_live_frame,
)
from fastmdxplora.gui.protein_preview import protein_preview_payload
from fastmdxplora.gui.server import (
    DashboardConfig,
    make_handler,
    start_dashboard_session,
    start_test_server,
)
from fastmdxplora.gui.structure_info import count_structure, ligand_atom_counts
from fastmdxplora.simulation.runner import _maybe_write_live_frame
from fastmdxplora.gui.telemetry import (
    TelemetryWriter,
    analyze_health,
    read_events,
    read_metrics,
    read_status,
)
from fastmdxplora.gui.trajectory_playback import (
    PlaybackUnavailable,
    neighborhood_residues,
    playback_info,
)


# ---------------------------------------------------------------------------
# Fixture: tiny four-atom helper used by several tests
# ---------------------------------------------------------------------------
def _tiny_pdb() -> str:
    return "\n".join(
        [
            "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 10.00           C",
            "ATOM      2  CA  GLY A   2       1.500   0.400   0.300  1.00 10.00           C",
            "ATOM      3  CA  SER A   3       2.800   1.200   0.000  1.00 10.00           C",
            "ATOM      4  CA  LEU A   4       4.100   1.400   0.900  1.00 10.00           C",
            "END",
        ]
    )


# ---------------------------------------------------------------------------
# Telemetry basics (unchanged behavior)
# ---------------------------------------------------------------------------
def test_telemetry_writer_creates_status_metrics_and_events(tmp_path: Path) -> None:
    simulation_dir = tmp_path / "run" / "simulation"
    writer = TelemetryWriter(simulation_dir, total_steps=100, planned_frames=10, timestep_fs=2.0, platform="CPU", target_temperature_K=300.0)
    writer.write_status(stage="production", status="running", current_step=25)
    writer.append_metric(stage="production", step=25, simulation_time_ns=0.00005, potential_energy=-123.4, total_energy=-100.0, temperature=299.0, progress_percent=25.0)
    writer.event("frame written")

    status = read_status(tmp_path / "run")
    metrics = read_metrics(tmp_path / "run")
    events = read_events(tmp_path / "run")

    assert status["stage"] == "production"
    assert status["current_step"] == 25
    assert metrics[-1]["potential_energy"] == "-123.4"
    assert events[-1]["message"] == "frame written"


def test_telemetry_writer_file_errors_do_not_raise(tmp_path: Path) -> None:
    blocked_path = tmp_path / "not_a_directory"
    blocked_path.write_text("occupied", encoding="utf-8")
    writer = TelemetryWriter(blocked_path)
    writer.write_status(stage="production")
    writer.append_metric(stage="production", step=1)
    writer.event("event")


def test_health_detects_nan_and_energy_spike() -> None:
    status = {"status": "running", "last_update_timestamp": "2999-01-01T00:00:00+00:00"}
    nan_health = analyze_health(status, [{"stage": "production", "total_energy": "nan"}])
    assert nan_health["state"] == "failed"
    assert "numerically unstable" in nan_health["explanation"]

    spike_health = analyze_health(status, [{"stage": "production", "total_energy": "-100.0"}, {"stage": "production", "total_energy": "50000.0"}])
    assert spike_health["state"] == "warning"
    assert "Energy increased sharply" in spike_health["explanation"]


def test_telemetry_readers_handle_missing_and_malformed_files(tmp_path: Path) -> None:
    run = tmp_path / "run"
    sim = run / "simulation"
    sim.mkdir(parents=True)

    assert read_status(run) == {}
    assert read_metrics(run) == []
    assert read_events(run) == []

    (sim / "live_status.json").write_text("{not-json", encoding="utf-8")
    (sim / "live_metrics.csv").write_text('timestamp,stage,total_energy\n"unterminated', encoding="utf-8")
    (sim / "live_events.log").write_text("free-form event without tabs\n", encoding="utf-8")

    assert read_status(run) == {}
    assert isinstance(read_metrics(run), list)
    assert read_events(run) == [{"timestamp": "", "level": "info", "message": "free-form event without tabs"}]


def test_analyze_health_temperature_and_stale_warning() -> None:
    stale = analyze_health({"status": "running", "last_update_timestamp": "2000-01-01T00:00:00+00:00"}, [], stale_after_seconds=1)
    assert stale["state"] == "warning"
    assert "stale" in stale["message"].lower()

    hot = analyze_health({"status": "running", "target_temperature_K": 300.0}, [{"temperature": "420.0"}])
    assert hot["state"] == "warning"
    assert "Temperature" in hot["message"]


# ---------------------------------------------------------------------------
# Protein preview behaviour (unchanged endpoint contract)
# ---------------------------------------------------------------------------
def test_protein_preview_unavailable_without_structure(tmp_path: Path) -> None:
    payload = protein_preview_payload(tmp_path / "run")
    assert payload["available"] is False
    assert payload["message"] == "No topology/PDB found yet."


def test_protein_preview_uses_existing_image_without_structure(tmp_path: Path) -> None:
    run = tmp_path / "run"
    preview = run / "report" / "dashboard_assets" / "protein_preview.png"
    preview.parent.mkdir(parents=True)
    preview.write_bytes(b"\x89PNG\r\n\x1a\n")

    payload = protein_preview_payload(run)
    assert payload["available"] is True
    assert payload["static_available"] is True
    assert payload["static_mode"] == "existing"
    assert payload["path"] == "report/dashboard_assets/protein_preview.png"


def test_structure_route_serves_topology(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setattr(protein_preview, "_find_pymol_executable", lambda _root: None)
    run = tmp_path / "run"
    topology = run / "simulation" / "topology.pdb"
    topology.parent.mkdir(parents=True)
    topology.write_text(_tiny_pdb(), encoding="utf-8")

    server, base_url = start_test_server(run)
    try:
        response = urlopen(f"{base_url}/structure/topology.pdb", timeout=5)
        body = response.read().decode("utf-8")
    finally:
        server.shutdown()
        server.server_close()

    assert response.headers["Cache-Control"] == "no-store"
    assert "ATOM" in body and "ALA" in body


def test_pymol_script_uses_padded_uncropped_camera(tmp_path: Path, monkeypatch) -> None:
    structure = tmp_path / "topology.pdb"
    output = tmp_path / "protein_preview.png"
    structure.write_text(_tiny_pdb(), encoding="utf-8")
    captured: dict[str, str] = {}

    def _fake_run(command, **_kwargs):
        script = Path(command[-1])
        captured["script"] = script.read_text(encoding="utf-8")
        output.write_bytes(b"\x89PNG\r\n\x1a\n")

    monkeypatch.setattr(protein_preview.subprocess, "run", _fake_run)
    protein_preview._render_with_pymol("/fake/pymol", structure, output)

    script = captured["script"]
    assert "viewport 1800, 1400" in script
    assert "show cartoon, prot" in script
    assert "spectrum count, rainbow, prot, byres=1" in script
    assert "set_color fastmdx_res_1" in script
    assert "color fastmdx_res_1, prot and resi 1 and chain A" in script
    assert "set ray_opaque_background, off" in script
    assert "set cartoon_fancy_helices, 1" in script
    assert "set cartoon_smooth_loops, 1" in script
    assert "set ray_trace_mode, 1" in script
    assert "set cartoon_tube_radius, 0.45" in script
    assert "set cartoon_sampling, 14" in script
    assert "orient prot" in script
    assert "center prot" in script
    assert "zoom prot, 1.8" in script
    assert "ray 1800, 1200" in script
    assert "width=1800" in script
    assert "height=1200" in script
    assert "dpi=300" in script


def test_gui_command_parses() -> None:
    args = _build_parser().parse_args(
        ["gui", "--output", "local_runs/my_run", "--port", "8766"]
    )
    assert args.command == "gui"
    assert args.output == "local_runs/my_run"
    assert args.port == 8766


# ---------------------------------------------------------------------------
# New dashboards: structure_info, ligand_detection, live_frames, trajectory_playback
# ---------------------------------------------------------------------------
def test_count_structure_basic_protein(tmp_path: Path) -> None:
    (tmp_path / "topology.pdb").write_text(_tiny_pdb(), encoding="utf-8")
    info = count_structure(tmp_path / "topology.pdb")
    assert info["valid"] is True
    assert info["n_chains"] == 1
    assert info["protein_residues"] == 4
    assert info["water_residues"] == 0
    assert info["ions"] == 0
    assert info["ligand_residues_total"] == 0


def test_count_structure_water_atom_record_not_double_counted(tmp_path: Path) -> None:
    # Many PDBs list waters as ATOM HOH. They must not appear in
    # protein_residues.
    pdb = "\n".join([
        "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N",
        "ATOM      2  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C",
        "ATOM      3  OH2 HOH B  10       5.000   5.000   5.000  1.00  0.00           O",
        "END",
    ])
    (tmp_path / "topology.pdb").write_text(pdb, encoding="utf-8")
    info = count_structure(tmp_path / "topology.pdb")
    assert info["valid"] is True
    assert info["protein_residues"] == 1
    assert info["water_residues"] == 1


def test_count_structure_detects_ligand(tmp_path: Path) -> None:
    pdb = "\n".join([
        "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N",
        "ATOM      2  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C",
        "HETATM    3  C1  EPE A 100       4.000   4.000   4.000  1.00 10.00           C",
        "HETATM    4  O1  EPE A 100       4.300   4.300   4.300  1.00 10.00           O",
        "END",
    ])
    (tmp_path / "topology.pdb").write_text(pdb, encoding="utf-8")
    info = count_structure(tmp_path / "topology.pdb")
    assert "EPE" in info["ligand_resnames"]


def test_count_structure_handles_missing(tmp_path: Path) -> None:
    info = count_structure(tmp_path / "no.pdb")
    assert info["valid"] is False


def test_count_structure_rejects_too_large_files(tmp_path: Path) -> None:
    big = tmp_path / "huge.pdb"
    big.write_text("HEADER\n", encoding="utf-8")
    big.write_bytes(b"X" * (9 * 1024 * 1024))  # > 8 MB
    info = count_structure(big)
    assert info["valid"] is False
    assert info["reason"] == "too-large"


def test_detect_ligands_basic(tmp_path: Path) -> None:
    pdb = "\n".join([
        "HETATM    1  C1  EPE A 100       4.000   4.000   4.000  1.00 10.00           C",
        "HETATM    2  C2  BEN B   1       5.000   5.000   5.000  1.00 10.00           C",
        "ATOM      3  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N",
        "END",
    ])
    (tmp_path / "topology.pdb").write_text(pdb, encoding="utf-8")
    info = count_structure(tmp_path / "topology.pdb")
    keys = [(ins["chain"], ins["resname"], ins["resi"]) for ins in info["ligand_instances"]]
    detected = detect_ligands(keys)
    assert set(detected["resnames"]) >= {"EPE"}


def test_detect_ligands_respects_explicit_override() -> None:
    keys = [("A", "NA", "1")]
    out = detect_ligands(keys, explicit=["NA"])
    assert out["resnames"] == ["NA"]
    assert out["instances"][0]["explicit"] is True


def test_detect_ligands_excludes_water_and_ions() -> None:
    keys = [("A", "HOH", "1"), ("A", "NA", "2"), ("A", "EPE", "3")]
    out = detect_ligands(keys)
    assert sorted(out["resnames"]) == ["EPE"]


def test_detect_ligands_includes_cofactors_when_requested() -> None:
    keys = [("A", "GOL", "1"), ("A", "EPE", "2")]
    default = detect_ligands(keys)
    with_cof = detect_ligands(keys, include_cofactors=True)
    assert "GOL" in with_cof["cofactors"]
    assert "GOL" not in default["resnames"]
    assert "GOL" in with_cof["resnames"]


def test_detect_ligands_lexicographic_sort(tmp_path: Path) -> None:
    # Single chain so the (chain, resname, resi) sort key collapses to
    # numeric resi ordering only.
    keys = [("A", "EPE", str(r)) for r in [10, 2, 100]]
    detected = detect_ligands(keys)
    assert [ins["resi"] for ins in detected["instances"]] == ["2", "10", "100"]


def test_ligand_atom_counts_distinguishes_ligands(tmp_path: Path) -> None:
    pdb = "\n".join([
        "HETATM    1  C1  EPE A 100       4.000   4.000   4.000  1.00 10.00           C",
        "HETATM    2  C2  EPE A 100       4.300   4.300   4.300  1.00 10.00           C",
        "HETATM    3  C1  BEN B   1       5.000   5.000   5.000  1.00 10.00           C",
        "HETATM    4  H1  HOH A   2       5.500   5.500   5.500  1.00  0.00           H",
        "END",
    ])
    (tmp_path / "topology.pdb").write_text(pdb, encoding="utf-8")
    counts = ligand_atom_counts(tmp_path / "topology.pdb")
    assert counts == {"EPE": 2, "BEN": 1}


def test_normalise_ligand_resname() -> None:
    assert normalise_ligand_resname("epe") == "EPE"
    assert normalise_ligand_resname("LIG ") == "LIG"
    assert normalise_ligand_resname("") is None
    assert normalise_ligand_resname(None) is None


def test_write_live_frame_round_trip(tmp_path: Path) -> None:
    sim = tmp_path / "simulation"
    result = write_live_frame(sim, pdb_text=_tiny_pdb(), frame_index=42)
    assert result["ok"] is True
    assert (sim / "live_frame.pdb").is_file()
    index = json.loads((sim / "live_frame_index.json").read_text(encoding="utf-8"))
    assert index["live_frame_available"] is True
    assert index["live_frame_index"] == 42


def test_write_live_frame_handles_missing_directory_safely(tmp_path: Path) -> None:
    # write_live_frame must not crash when the target "directory" is a
    # regular file. Putting a regular file at the target slot makes
    # Path.mkdir(parents=True, exist_ok=True) fail both on POSIX and
    # Windows. Cross-platform safe; the original /nonexistent/... path
    # was environment-dependent on Windows because C:\ is normally
    # writable so the parent.mkdir succeeded.
    blocked = tmp_path / "i_am_a_file"
    blocked.write_text("not a directory", encoding="utf-8")
    result = write_live_frame(blocked, pdb_text="ATOM\n", frame_index=1)
    assert result["ok"] is False


def test_playback_returns_unavailable_when_missing(tmp_path: Path) -> None:
    info = playback_info(tmp_path)
    assert info["playback_available"] is False


def test_playback_generation_failure_is_safe(tmp_path: Path, monkeypatch) -> None:
    # Force a failure inside the writer.
    import fastmdxplora.gui.trajectory_playback as mod
    monkeypatch.setattr(mod, "_import_mdtraj", lambda: None)
    (tmp_path / "simulation").mkdir(parents=True)
    (tmp_path / "simulation" / "topology.pdb").write_text(_tiny_pdb(), encoding="utf-8")
    (tmp_path / "simulation" / "production.dcd").write_bytes(b"\x00\x00")
    info = playback_info(tmp_path)
    assert info["playback_available"] is False
    assert info["reason"] == "mdtraj-not-installed"


def test_neighborhood_residues_per_atom_check(tmp_path: Path) -> None:
    # Atom-wide coordinates in standard PDB column widths so the parser
    # finds the values in the right slots (cols 31-38, 39-46, 47-54).
    pdb = "\n".join([
        "HETATM    1  C1  EPE A 100       4.000   4.000   4.000  1.00 10.00           C",
        "ATOM      2  N   ALA A   1       4.500   4.500   4.500  1.00  0.00           N",  # near
        "ATOM      3  N   ALA A   2      30.000  30.000  30.000  1.00  0.00           N",  # far
        "END",
    ])
    (tmp_path / "topology.pdb").write_text(pdb, encoding="utf-8")
    residues = neighborhood_residues(topology_path=tmp_path / "topology.pdb", ligand_resname="EPE", cutoff_angstrom=5.0)
    assert ("A", 1) in residues
    assert ("A", 2) not in residues



def test_telemetry_stage_states_survive_multiple_writers(tmp_path: Path) -> None:
    sim = tmp_path / "run" / "simulation"
    orchestrator_writer = TelemetryWriter(sim)
    runner_writer = TelemetryWriter(sim)
    orchestrator_writer.write_status(
        stage="setup",
        status="running",
        stage_states={
            "setup": "current", "minimization": "waiting", "nvt": "waiting",
            "npt": "waiting", "production": "waiting", "analysis": "waiting",
            "report": "waiting",
        },
    )
    runner_writer.mark_stage("minimization", "current", status="running")
    runner_writer.mark_stage("minimization", "completed", status="running")
    orchestrator_writer.mark_stage("analysis", "current", status="running")
    status = read_status(tmp_path / "run")
    assert status["stage_states"]["setup"] == "current"
    assert status["stage_states"]["minimization"] == "completed"
    assert status["stage_states"]["analysis"] == "current"


def test_live_frame_history_builds_playback(tmp_path: Path) -> None:
    sim = tmp_path / "simulation"
    write_live_frame(
        sim, pdb_text=_tiny_pdb(), frame_index=10, stage="nvt",
        simulation_time_ns=0.001, archive=True,
    )
    moved = _tiny_pdb().replace("   0.000   0.000   0.000", "   0.200   0.000   0.000")
    write_live_frame(
        sim, pdb_text=moved, frame_index=20, stage="nvt",
        simulation_time_ns=0.002, archive=True,
    )
    history = read_live_frame_history(sim)
    assert history["count"] == 2
    info = playback_info(tmp_path, max_browser_frames=20)
    assert info["playback_available"] is True
    assert info["source_kind"] == "live-history"
    assert info["n_frames_browser"] == 2
    text = (sim / "playback.pdb").read_text(encoding="utf-8")
    assert text.count("MODEL") == 2


def test_dashboard_display_pdb_strips_solvent_and_keeps_ligand() -> None:
    pdb = "\n".join([
        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 10.00           C",
        "HETATM    2  O   HOH B   2       2.000   2.000   2.000  1.00 10.00           O",
        "HETATM    3  NA   NA C   3       3.000   3.000   3.000  1.00 10.00          NA",
        "HETATM    4  C1  LIG D   4       4.000   4.000   4.000  1.00 10.00           C",
        "END",
    ])
    display = dashboard_display_pdb(pdb)
    assert "ALA" in display
    assert "LIG" in display
    assert "HOH" not in display
    assert " NA C" not in display

# ---------------------------------------------------------------------------
# Endpoint contract: the redesigned dashboard's HTML/JS layer
# ---------------------------------------------------------------------------
def test_dashboard_html_uses_separated_template() -> None:
    """The HTML shell lives on disk next to the module rather than as a giant f-string."""
    template_path = Path(make_handler.__module__.split(".")[0])
    import fastmdxplora.gui as gui_pkg
    layout = Path(gui_pkg.__file__).with_name("templates") / "dashboard.html"
    assert layout.is_file()


def test_dashboard_html_has_aai_branding(tmp_path: Path) -> None:
    run = tmp_path / "run"
    writer = TelemetryWriter(run / "simulation")
    writer.write_status(stage="NVT")
    server, base_url = start_test_server(run)
    try:
        html = urlopen(f"{base_url}/", timeout=5).read().decode("utf-8")
    finally:
        server.shutdown()
        server.server_close()

    # Branding hierarchy
    assert "FastMDXplora" in html
    assert "Fully Automated SysTem for Molecular Dynamics eXploration" in html

    # Sidebar / top bar / nav structure
    assert "sidebar" in html
    assert "top-bar" in html
    assert "data-view-link" in html

    # Asset wiring
    assert "/static/dashboard.css" in html
    assert "/static/dashboard.js" in html
    assert "/static/charts.js" in html
    assert "/static/molecule-viewer.js" in html
    assert "/static/3Dmol-min.js" in html

    # Pages
    for page in ("overview", "live", "viewer", "analysis", "files", "settings"):
        assert f'data-page="{page}"' in html

    # Loading screen + branded particles
    assert "loading-screen" in html

    # No CDN references
    assert "cdn" not in html.lower()
    assert "googleapis" not in html.lower()
    # The documentation + GitHub sidebar links are the only legitimate
    # external refs, both harmless but explicit.
    for allowed in (
        "https://fastmdxplora.readthedocs.io",
        "https://github.com/aai-research-lab/FastMDXplora",
    ):
        assert allowed in html
    assert html.count("https://") == 2


def test_dashboard_css_uses_black_scientific_palette() -> None:
    import fastmdxplora.gui as gui_pkg

    static = Path(gui_pkg.__file__).with_name("static")
    # Design tokens live in theme.css (shared with the static report);
    # dashboard.css carries the layout that consumes them.
    theme = (static / "theme.css").read_text(encoding="utf-8")
    css = (static / "dashboard.css").read_text(encoding="utf-8")
    for token in ("--background-primary", "--accent-cyan", "--accent-violet",
                  "--text-primary", "Inter"):
        assert token in theme
    assert "prefers-reduced-motion" in css


def test_dashboard_endpoint_no_inline_html_css() -> None:
    """The server.py module no longer contains the giant HTML f-string.

    We assert the absence of known-fingerprinted fragments from the old
    inline template — a tautology-free check.
    """
    server_path = Path(make_handler.__module__.split(".")[0])
    import fastmdxplora.gui as gui_pkg
    src = (Path(gui_pkg.__file__).with_name("server.py")).read_text(encoding="utf-8")
    for marker in (
        'onclick="setProteinPreviewExpanded',
        'spectrum count, rainbow, prot, byres=1',
        'function setProteinPreviewExpanded(expanded)',
    ):
        assert marker not in src, f"legacy inline marker found: {marker!r}"


def test_dashboard_server_serves_new_routes(tmp_path: Path) -> None:
    run = tmp_path / "run"
    sim = run / "simulation"
    sim.mkdir(parents=True)
    (sim / "topology.pdb").write_text(_tiny_pdb(), encoding="utf-8")
    write_live_frame(sim, pdb_text=_tiny_pdb(), frame_index=10)

    server, base_url = start_test_server(run)
    try:
        for path in (
            "/api/status",
            "/api/metrics",
            "/api/events",
            "/api/artifacts",
            "/api/files",          # alias
            "/api/results",
            "/api/analyses",       # alias
            "/api/protein-preview",
            "/api/structure-info",
            "/api/ligands",
            "/api/live-frame-index",
            "/api/live-coordinates",
            "/api/playback-info",
            "/structure/topology.pdb",
            "/structure/live-frame.pdb",
        ):
            response = urlopen(f"{base_url}{path}", timeout=5)
            response.read()
            assert response.headers["Cache-Control"] == "no-store", path
    finally:
        server.shutdown()
        server.server_close()


def test_structure_info_endpoint_payload_keys(tmp_path: Path) -> None:
    run = tmp_path / "run"
    sim = run / "simulation"
    sim.mkdir(parents=True)
    (sim / "topology.pdb").write_text(_tiny_pdb(), encoding="utf-8")
    server, base_url = start_test_server(run, config=DashboardConfig(ligand_resname="LIG"))
    try:
        payload = json.loads(urlopen(f"{base_url}/api/structure-info", timeout=5).read())
    finally:
        server.shutdown()
        server.server_close()
    assert payload["valid"] is True
    assert payload["n_chains"] == 1
    assert payload["protein_residues"] == 4
    assert payload["explicit_ligand"] == "LIG"


def test_live_frame_endpoint_serves_pdb(tmp_path: Path) -> None:
    run = tmp_path / "run"
    sim = run / "simulation"
    sim.mkdir(parents=True)
    write_live_frame(sim, pdb_text=_tiny_pdb(), frame_index=5)
    server, base_url = start_test_server(run)
    try:
        body = urlopen(f"{base_url}/structure/live-frame.pdb", timeout=5).read().decode("utf-8")
    finally:
        server.shutdown()
        server.server_close()
    assert "ATOM" in body and "ALA" in body


def test_live_frame_endpoint_404_when_missing(tmp_path: Path) -> None:
    run = tmp_path / "run"
    sim = run / "simulation"
    sim.mkdir(parents=True)
    server, base_url = start_test_server(run)
    try:
        try:
            urlopen(f"{base_url}/structure/live-frame.pdb", timeout=5)
        except HTTPError as exc:
            assert exc.code == 404
        else:
            raise AssertionError("expected 404 on missing live frame")
    finally:
        server.shutdown()
        server.server_close()


def test_playback_info_unavailable_when_no_trajectory(tmp_path: Path) -> None:
    run = tmp_path / "run"
    run.mkdir()
    server, base_url = start_test_server(run)
    try:
        payload = json.loads(urlopen(f"{base_url}/api/playback-info", timeout=5).read())
    finally:
        server.shutdown()
        server.server_close()
    assert payload["playback_available"] is False


def test_artifacts_response_includes_files(tmp_path: Path) -> None:
    run = tmp_path / "run"
    sim = run / "simulation"
    sim.mkdir(parents=True)
    (sim / "live_status.json").write_text("{}", encoding="utf-8")
    server, base_url = start_test_server(run)
    try:
        payload = json.loads(urlopen(f"{base_url}/api/files", timeout=5).read())
    finally:
        server.shutdown()
        server.server_close()
    assert any(record["path"].endswith("live_status.json") for record in payload["artifacts"])


def test_results_response_keys_match_dashboard_bindings(tmp_path: Path) -> None:
    run = tmp_path / "run"
    sim = run / "simulation"
    sim.mkdir(parents=True)
    (sim / "live_status.json").write_text("{}", encoding="utf-8")
    server, base_url = start_test_server(run)
    try:
        payload = json.loads(urlopen(f"{base_url}/api/results", timeout=5).read())
    finally:
        server.shutdown()
        server.server_close()
    for key in ("refreshed_at", "has_analysis", "has_report", "summary", "system", "phases", "plots", "key_plots", "reports", "artifacts"):
        assert key in payload


def test_dashboard_session_uses_next_port_when_busy(tmp_path: Path) -> None:
    run = tmp_path / "run"
    run.mkdir()
    blocker = socket.socket()
    blocker.bind(("127.0.0.1", 0))
    blocker.listen()
    requested_port = blocker.getsockname()[1]
    session = start_dashboard_session(output=run, host="127.0.0.1", port=requested_port, max_port_tries=5)
    try:
        assert session.requested_port == requested_port
        assert session.port != requested_port
        assert session.port_was_changed is True
    finally:
        session.stop()
        blocker.close()


def test_dashboard_session_stop_is_safe(tmp_path: Path) -> None:
    run = tmp_path / "run"
    run.mkdir()
    session = start_dashboard_session(output=run, port=0)
    assert session.thread.is_alive()
    session.stop()
    assert not session.thread.is_alive()


def test_live_server_rejects_artifact_path_traversal(tmp_path: Path) -> None:
    run = tmp_path / "run"
    run.mkdir()
    outside = tmp_path / "outside.txt"
    outside.write_text("secret", encoding="utf-8")
    server, base_url = start_test_server(run)
    try:
        try:
            urlopen(f"{base_url}/artifacts/../outside.txt", timeout=5)
        except HTTPError as exc:
            assert exc.code in {403, 404}
        else:
            raise AssertionError("path traversal unexpectedly succeeded")
    finally:
        server.shutdown()
        server.server_close()


def test_live_server_does_not_serve_static_path_traversal(tmp_path: Path) -> None:
    run = tmp_path / "run"
    run.mkdir()
    server, base_url = start_test_server(run)
    try:
        try:
            urlopen(f"{base_url}/static/../server.py", timeout=5)
        except HTTPError as exc:
            assert exc.code == 404
        else:
            raise AssertionError("static traversal unexpectedly succeeded")
    finally:
        server.shutdown()
        server.server_close()


def test_static_3dmol_asset_is_served_locally(tmp_path: Path) -> None:
    run = tmp_path / "run"
    run.mkdir()
    server, base_url = start_test_server(run)
    try:
        response = urlopen(f"{base_url}/static/3Dmol-min.js", timeout=5)
        body = response.read().decode("utf-8", errors="ignore")
    finally:
        server.shutdown()
        server.server_close()
    assert response.headers["Cache-Control"] == "no-store"
    assert "$3Dmol" in body
    assert "http://" not in body[:1000]
    assert "https://" not in body[:1000]


def test_live_json_endpoints_are_no_store(tmp_path: Path) -> None:
    run = tmp_path / "run"
    writer = TelemetryWriter(run / "simulation")
    writer.write_status(stage="production")
    server, base_url = start_test_server(run)
    try:
        for endpoint in (
            "status", "metrics", "events", "artifacts", "results", "files",
            "protein-preview", "structure-info", "ligands",
            "live-frame-index", "live-coordinates", "playback-info",
            "analyses",
        ):
            response = urlopen(f"{base_url}/api/{endpoint}", timeout=5)
            response.read()
            assert response.headers["Cache-Control"] == "no-store", endpoint
    finally:
        server.shutdown()
        server.server_close()


def test_live_json_endpoints_are_sane_for_empty_run(tmp_path: Path) -> None:
    run = tmp_path / "empty_run"
    server, base_url = start_test_server(run)
    try:
        status_payload = json.loads(urlopen(f"{base_url}/api/status", timeout=5).read())
        metrics_payload = json.loads(urlopen(f"{base_url}/api/metrics", timeout=5).read())
        events_payload = json.loads(urlopen(f"{base_url}/api/events", timeout=5).read())
        artifacts_payload = json.loads(urlopen(f"{base_url}/api/artifacts", timeout=5).read())
        results_payload = json.loads(urlopen(f"{base_url}/api/results", timeout=5).read())
        structure_payload = json.loads(urlopen(f"{base_url}/api/structure-info", timeout=5).read())
        ligands_payload = json.loads(urlopen(f"{base_url}/api/ligands", timeout=5).read())
    finally:
        server.shutdown()
        server.server_close()
    assert status_payload["status"] == {}
    assert status_payload["health"]["state"] == "unknown"
    assert metrics_payload["metrics"] == []
    assert events_payload["events"] == []
    assert artifacts_payload["artifacts"] == []
    assert results_payload["has_analysis"] is False
    assert results_payload["has_report"] is False
    assert results_payload["plots"] == []
    assert structure_payload["valid"] is False
    assert ligands_payload["valid"] is False


def test_static_assets_are_packaged(tmp_path: Path) -> None:
    import fastmdxplora.gui as gui_pkg
    static = Path(gui_pkg.__file__).with_name("static")
    assert (static / "dashboard.css").is_file()
    assert (static / "dashboard.js").is_file()
    assert (static / "charts.js").is_file()
    assert (static / "molecule-viewer.js").is_file()


def test_static_endpoint_serves_packaged_assets(tmp_path: Path) -> None:
    run = tmp_path / "run"
    run.mkdir()
    server, base_url = start_test_server(run)
    try:
        response = urlopen(f"{base_url}/static/dashboard.css", timeout=5)
        body = response.read().decode("utf-8", errors="ignore")
    finally:
        server.shutdown()
        server.server_close()
    assert response.headers["Content-Type"].startswith("text/css")
    assert ".sidebar" in body

# ---------------------------------------------------------------------------
# Regression coverage for the functional dashboard wiring
# ---------------------------------------------------------------------------
def test_find_structure_prefers_prepared_solute_over_solvated_topology(
    tmp_path: Path,
) -> None:
    run = tmp_path / "run"
    prepared = run / "setup" / "prepared.pdb"
    topology = run / "simulation" / "topology.pdb"
    solvated = run / "setup" / "solvated.pdb"
    for path, marker in (
        (prepared, "PREPARED"),
        (topology, "TOPOLOGY"),
        (solvated, "SOLVATED"),
    ):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(f"REMARK {marker}\n{_tiny_pdb()}\n", encoding="utf-8")

    assert protein_preview.find_structure(run) == prepared


def test_read_metrics_falls_back_to_openmm_energy_csv(tmp_path: Path) -> None:
    sim = tmp_path / "run" / "simulation"
    sim.mkdir(parents=True)
    (sim / "energy.csv").write_text(
        '#"Step","Time (ps)","Potential Energy (kJ/mole)",'
        '"Kinetic Energy (kJ/mole)","Total Energy (kJ/mole)",'
        '"Temperature (K)","Box Volume (nm^3)","Density (g/mL)",'
        '"Speed (ns/day)","Progress (%)"\n'
        '100,2.0,-1000.5,200.5,-800.0,299.5,12.0,0.997,8.25,50.0\n',
        encoding="utf-8",
    )

    metrics = read_metrics(tmp_path / "run")

    assert len(metrics) == 1
    assert metrics[0]["step"] == "100"
    assert metrics[0]["simulation_time_ns"] == "0.002"
    assert metrics[0]["potential_energy"] == "-1000.5"
    assert metrics[0]["temperature"] == "299.5"
    assert metrics[0]["density"] == "0.997"
    assert metrics[0]["speed"] == "8.25"


def test_ligands_endpoint_returns_detected_instances(tmp_path: Path) -> None:
    run = tmp_path / "run"
    prepared = run / "setup" / "prepared.pdb"
    prepared.parent.mkdir(parents=True)
    prepared.write_text(
        "\n".join(
            [
                "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N",
                "HETATM    2  C1  EPE B 101       4.000   4.000   4.000  1.00 10.00           C",
                "HETATM    3  O1  EPE B 101       4.300   4.300   4.300  1.00 10.00           O",
                "END",
            ]
        ),
        encoding="utf-8",
    )

    server, base_url = start_test_server(run)
    try:
        payload = json.loads(urlopen(f"{base_url}/api/ligands", timeout=5).read())
    finally:
        server.shutdown()
        server.server_close()

    assert payload["valid"] is True
    assert payload["resnames"] == ["EPE"]
    assert payload["ligands"] == [
        {"chain": "B", "resname": "EPE", "resi": "101", "explicit": False}
    ]


def test_results_endpoint_discovers_analysis_and_report_outputs(
    tmp_path: Path,
) -> None:
    run = tmp_path / "run"
    analysis = run / "analysis" / "rmsd"
    report = run / "report"
    analysis.mkdir(parents=True)
    report.mkdir(parents=True)

    (run / "analysis" / "analysis_manifest.json").write_text(
        json.dumps(
            {
                "n_frames": 10,
                "results": {
                    "rmsd": {
                        "status": "ok",
                        "message": "RMSD complete",
                    },
                    "rg": {
                        "status": "skipped",
                        "message": "Not enough frames",
                    },
                },
            }
        ),
        encoding="utf-8",
    )
    (analysis / "rmsd.png").write_bytes(b"\x89PNG\r\n\x1a\n")
    (analysis / "rmsd.csv").write_text("frame,rmsd\n0,0.0\n", encoding="utf-8")
    (report / "dashboard.html").write_text("<html></html>", encoding="utf-8")
    (report / "report.md").write_text("# Report\n", encoding="utf-8")
    (report / "slides.pptx").write_bytes(b"PK")
    (report / "project_bundle.zip").write_bytes(b"PK")

    server, base_url = start_test_server(run)
    try:
        payload = json.loads(urlopen(f"{base_url}/api/results", timeout=5).read())
    finally:
        server.shutdown()
        server.server_close()

    assert payload["has_analysis"] is True
    assert payload["has_report"] is True
    assert {item["name"] for item in payload["analyses"]} == {"rmsd", "rg"}
    rmsd = next(item for item in payload["analyses"] if item["name"] == "rmsd")
    assert rmsd["plot"]["path"] == "analysis/rmsd/rmsd.png"
    assert any(item["path"] == "analysis/rmsd/rmsd.csv" for item in rmsd["artifacts"])
    assert [item["path"] for item in payload["reports"]] == [
        "report/dashboard.html",
        "report/report.md",
        "report/slides.pptx",
        "report/project_bundle.zip",
    ]


def test_artifact_download_route_sets_attachment_header(tmp_path: Path) -> None:
    run = tmp_path / "run"
    report = run / "report"
    report.mkdir(parents=True)
    (report / "report.md").write_text("# Report\n", encoding="utf-8")

    server, base_url = start_test_server(run)
    try:
        response = urlopen(
            f"{base_url}/artifacts/report/report.md?download=1",
            timeout=5,
        )
        response.read()
    finally:
        server.shutdown()
        server.server_close()

    assert response.headers["Content-Disposition"] == 'attachment; filename="report.md"'


def test_frontend_assets_wire_functional_dashboard_sections() -> None:
    import fastmdxplora.gui as gui_pkg

    root = Path(gui_pkg.__file__).parent
    html = (root / "templates" / "dashboard.html").read_text(encoding="utf-8")
    css = (root / "static" / "dashboard.css").read_text(encoding="utf-8")
    dashboard_js = (root / "static" / "dashboard.js").read_text(encoding="utf-8")
    charts_js = (root / "static" / "charts.js").read_text(encoding="utf-8")
    viewer_js = (root / "static" / "molecule-viewer.js").read_text(encoding="utf-8")

    assert 'id="reports-files"' in html
    assert 'id="simulation-files"' in html
    assert 'id="analysis-files"' in html
    assert 'id="mini-preview-canvas"' in html
    assert "[hidden]" in css and "display: none !important" in css
    assert "/api/results" in dashboard_js
    assert "renderFiles" in dashboard_js
    assert "renderAnalysis" in dashboard_js
    assert "ResizeObserver" in charts_js
    assert "mini-preview-canvas" in viewer_js
    assert "/structure/topology.pdb" in viewer_js


def test_dashboard_refresh_seconds_are_injected_into_html(tmp_path: Path) -> None:
    run = tmp_path / "run"
    run.mkdir()
    server, base_url = start_test_server(
        run,
        config=DashboardConfig(refresh_seconds=1.5),
    )
    try:
        html = urlopen(base_url, timeout=5).read().decode("utf-8")
    finally:
        server.shutdown()
        server.server_close()
    assert 'data-refresh-seconds="1.5"' in html
    assert 'id="setting-refresh-seconds" min="1" max="60" step="1" value="1.5"' in html


def test_runner_live_frame_helper_calls_writer_with_valid_keyword(
    tmp_path: Path, monkeypatch
) -> None:
    calls: list[dict] = []

    def fake_write_openmm_live_frame(output_dir, **kwargs):
        calls.append({"output_dir": output_dir, **kwargs})
        return {"ok": True}

    import fastmdxplora.gui.live_frames as live_frames_module

    monkeypatch.setattr(
        live_frames_module,
        "write_openmm_live_frame",
        fake_write_openmm_live_frame,
    )

    class FakeState:
        def getPositions(self):
            return "positions"

    class FakeContext:
        def getState(self, **_kwargs):
            return FakeState()

    simulation = SimpleNamespace(context=FakeContext(), topology="topology")
    telemetry = SimpleNamespace(root=tmp_path / "simulation")
    omm = {"PDBFile": SimpleNamespace(writeFile=lambda *_args, **_kwargs: None)}

    _maybe_write_live_frame(
        omm,
        simulation,
        telemetry,
        250,
        stage="nvt",
        simulation_time_ns=0.0005,
    )

    assert len(calls) == 1
    assert calls[0]["pdbfile_writer"] is omm["PDBFile"].writeFile
    assert calls[0]["frame_index"] == 250
    assert calls[0]["stage"] == "nvt"
    assert calls[0]["archive"] is True


def test_dashboard_frame_interval_is_forwarded_to_explore_config() -> None:
    config = {"simulation": {"temperature_K": 300.0}}
    args = SimpleNamespace(dashboard_frame_interval=125)
    _enable_dashboard_telemetry(config, args)
    assert config["simulation"]["live_telemetry"] is True
    assert config["simulation"]["telemetry_interval"] == 125
    assert config["simulation"]["temperature_K"] == 300.0


def test_results_endpoint_pairs_png_preview_with_svg_download(tmp_path: Path) -> None:
    run = tmp_path / "run"
    analysis = run / "analysis" / "rg"
    analysis.mkdir(parents=True)
    (run / "analysis" / "analysis_manifest.json").write_text(
        json.dumps({"results": {"rg": {"status": "ok", "message": "rg: ok"}}}),
        encoding="utf-8",
    )
    (analysis / "rg.png").write_bytes(b"\x89PNG\r\n\x1a\n")
    (analysis / "rg.svg").write_text(
        '<svg xmlns="http://www.w3.org/2000/svg"></svg>', encoding="utf-8"
    )
    (analysis / "rg.dat").write_text("frame,rg_nm\n0,1.0\n", encoding="utf-8")

    server, base_url = start_test_server(run)
    try:
        payload = json.loads(urlopen(f"{base_url}/api/results", timeout=5).read())
    finally:
        server.shutdown()
        server.server_close()

    record = payload["analyses"][0]
    assert record["plot"]["path"] == "analysis/rg/rg.png"
    assert record["plot"]["svg_path"] == "analysis/rg/rg.svg"
    assert record["plot"]["svg_download_href"].startswith("/artifacts/analysis/rg/rg.svg")
    assert payload["svg_figure_count"] == 1


def test_svg_bundle_endpoint_downloads_generated_vector_figures(tmp_path: Path) -> None:
    import io
    import zipfile

    run = tmp_path / "run"
    figure = run / "analysis" / "rmsd" / "rmsd.svg"
    figure.parent.mkdir(parents=True)
    figure.write_text('<svg xmlns="http://www.w3.org/2000/svg"></svg>', encoding="utf-8")

    server, base_url = start_test_server(run)
    try:
        response = urlopen(f"{base_url}/analysis-figures-svg.zip", timeout=5)
        data = response.read()
    finally:
        server.shutdown()
        server.server_close()

    assert response.headers["Content-Type"] == "application/zip"
    assert "fastmdxplora_svg_figures.zip" in response.headers["Content-Disposition"]
    with zipfile.ZipFile(io.BytesIO(data)) as archive:
        assert archive.namelist() == ["analysis/rmsd/rmsd.svg"]


def test_dashboard_assets_include_svg_download_and_first_model_miniviewer_fix() -> None:
    root = Path(__file__).resolve().parents[1]
    html = (root / "src" / "fastmdxplora" / "gui" / "templates" / "dashboard.html").read_text(encoding="utf-8")
    dashboard_js = (root / "src" / "fastmdxplora" / "gui" / "static" / "dashboard.js").read_text(encoding="utf-8")
    viewer_js = (root / "src" / "fastmdxplora" / "gui" / "static" / "molecule-viewer.js").read_text(encoding="utf-8")

    assert 'id="download-all-svg"' in html
    assert "Download SVG" in dashboard_js
    assert "const hadModel" in viewer_js
    assert "opts.center !== false || !hadModel" in viewer_js
    assert "resolveProteinSelection" in viewer_js


def test_run_dependent_nav_sections_are_marked(tmp_path: Path) -> None:
    """Sections that need a run are marked so the shell can dim them.

    Opening the GUI in an empty workspace should point at the builder rather
    than offering five views that have nothing to show.
    """
    run = tmp_path / "run"
    run.mkdir()
    server, base_url = start_test_server(run)
    try:
        html = urlopen(f"{base_url}/", timeout=5).read().decode("utf-8")
    finally:
        server.shutdown()
        server.server_close()

    for view in ("overview", "live", "viewer", "analysis", "files"):
        assert f'data-view-link="{view}" data-requires-run="true"' in html
    # The builder is always reachable; it is what an empty workspace needs.
    assert 'data-view-link="builder" data-requires-run' not in html


def test_live_gui_and_static_report_share_one_theme(tmp_path: Path) -> None:
    """Both surfaces read their design tokens from ``live/static/theme.css``.

    The live GUI links the file; a generated report inlines it so the page
    stays self-contained. Duplicating the palette is what made the two look
    like different products, so this asserts a single source of truth.
    """
    import datetime

    import fastmdxplora.gui as gui_pkg
    from fastmdxplora.gui import report_dashboard

    theme = Path(gui_pkg.__file__).with_name("static") / "theme.css"
    assert theme.is_file(), "theme.css must ship with the package"
    tokens = theme.read_text(encoding="utf-8")
    assert "--background-primary" in tokens
    assert "--accent-cyan" in tokens

    # The GUI links it ahead of its own stylesheet.
    template = Path(gui_pkg.__file__).with_name("templates") / "dashboard.html"
    shell = template.read_text(encoding="utf-8")
    assert shell.index("/static/theme.css") < shell.index("/static/dashboard.css")

    # The static report inlines the very same tokens.
    html = report_dashboard._render_dashboard(
        title="Theme check",
        system="1L2Y",
        status="completed",
        generated_at=datetime.datetime.now(datetime.timezone.utc),
        phase_notice="",
        cards=[],
        sections=[],
        links=[],
        phase_rows=[],
        metrics=[],
        output_folder=str(tmp_path),
        live_html="",
    )
    assert "--background-primary" in html
    assert "--accent-cyan" in html
    # The superseded navy palette must not reappear.
    assert "#07101b" not in html
    assert "var(--bg)" not in html


def test_results_payload_carries_report_panels(tmp_path: Path) -> None:
    """The GUI and the static report show the same numbers.

    Summary cards, metric statistics, and phase rows are computed once in
    ``gui.report_dashboard`` and served to the browser, rather than the two
    surfaces calculating them independently and drifting.
    """
    import json

    from fastmdxplora.gui.server import _results_payload

    run = tmp_path / "run"
    (run / "analysis").mkdir(parents=True)
    (run / "simulation").mkdir(parents=True)
    (run / "manifest.json").write_text(
        json.dumps(
            {
                "system": "1L2Y",
                "phases": [
                    {
                        "name": "setup",
                        "status": "ok",
                        "started_at": "2026-08-01T00:00:00+00:00",
                        "finished_at": "2026-08-01T00:00:06+00:00",
                    },
                    {
                        "name": "simulation",
                        "status": "ok",
                        "started_at": "2026-08-01T00:00:06+00:00",
                        "finished_at": "2026-08-01T00:01:06+00:00",
                    },
                ],
            }
        ),
        encoding="utf-8",
    )
    (run / "analysis" / "analysis_manifest.json").write_text(
        json.dumps({"n_frames": 8, "n_atoms": 5080}), encoding="utf-8"
    )
    (run / "simulation" / "simulation_parameters.json").write_text(
        json.dumps({"parameters": {"temperature_K": 300.0, "production_steps": 400}}),
        encoding="utf-8",
    )

    payload = _results_payload(run)

    labels = {card["label"] for card in payload["summary_cards"]}
    assert {"Project status", "Phases completed", "Frames", "Atom count"} <= labels

    phases = {row["name"]: row["status"] for row in payload["phase_rows"]}
    assert phases["Setup"] == "ok"
    assert phases["Analysis"] == "not-run"

    metrics = {row["metric"] for row in payload["metric_rows"]}
    assert "Frame count" in metrics


def test_report_panels_never_break_the_dashboard(tmp_path: Path) -> None:
    """A malformed run must not take the whole payload down."""
    from fastmdxplora.gui.server import _results_payload

    run = tmp_path / "broken"
    run.mkdir()
    (run / "manifest.json").write_text("{not json", encoding="utf-8")

    payload = _results_payload(run)
    assert payload["summary_cards"] == [] or isinstance(payload["summary_cards"], list)
    assert "artifacts" in payload


def test_overview_hosts_the_report_panels(tmp_path: Path) -> None:
    """Overview has containers for the summary cards, phases, and statistics.

    These mirror what the generated report shows, so the browser is no longer
    the poorer of the two surfaces.
    """
    run = tmp_path / "run"
    run.mkdir()
    server, base_url = start_test_server(run)
    try:
        html = urlopen(f"{base_url}/", timeout=5).read().decode("utf-8")
    finally:
        server.shutdown()
        server.server_close()

    for container in (
        'id="overview-summary-cards"',
        'id="overview-phase-rows"',
        'id="overview-stat-rows"',
    ):
        assert container in html

    # Both tables carry the statistics columns the report shows.
    assert "Std. dev." in html
    assert "Trajectory statistics" in html
    assert "Exploration progress" in html

    import fastmdxplora.gui as gui_pkg

    script = (Path(gui_pkg.__file__).with_name("static") / "dashboard.js").read_text(
        encoding="utf-8"
    )
    assert "function renderReportPanels(" in script
    assert "renderReportPanels(payload);" in script


def test_results_payload_carries_categorised_sections(tmp_path: Path) -> None:
    """The GUI receives the report's grouping, not just a flat figure list."""
    import json

    from fastmdxplora.gui.server import _results_payload

    run = tmp_path / "run"
    for folder in ("analysis/rmsd", "analysis/cluster", "analysis/dimred", "report"):
        (run / folder).mkdir(parents=True)
    (run / "manifest.json").write_text(json.dumps({"system": "1L2Y", "phases": []}))
    (run / "analysis" / "analysis_manifest.json").write_text(json.dumps({"n_frames": 8}))
    for rel in (
        "analysis/rmsd/rmsd.png",
        "analysis/cluster/kmeans_trajectory.png",
        "analysis/dimred/pca.png",
    ):
        (run / rel).write_bytes(b"\x89PNG\r\n\x1a\n")

    payload = _results_payload(run)
    sections = payload["analysis_sections"]
    assert sections, "expected categorised sections"

    titles = {section["title"] for section in sections}
    # Each analysis keeps its own heading rather than being pooled.
    assert {"RMSD", "Clustering", "Dimensionality Reduction"} <= titles

    panels = [panel for section in sections for panel in section["panels"]]
    assert any(panel["title"] == "RMSD" for panel in panels)
    # Hrefs must be servable by the GUI, not relative to a report folder.
    assert all(panel["href"].startswith("/artifacts/") for panel in panels)


def test_data_files_are_summarised_without_plotting(tmp_path: Path) -> None:
    """Panel captions come straight from the data now.

    They used to be produced as a side effect of rendering a second copy of
    every figure, which meant drawing 17 charts nothing displayed.
    """
    from fastmdxplora.gui.report_dashboard import _summarise_data_file

    series = tmp_path / "rmsd.dat"
    series.write_text("\n".join(f"{i} {0.2 + i * 0.01:.4f}" for i in range(10)),
                      encoding="utf-8")
    assert _summarise_data_file(series, "line").startswith("avg ")

    clusters = tmp_path / "clusters.dat"
    clusters.write_text("\n".join(f"{i} {i % 3}" for i in range(30)), encoding="utf-8")
    assert _summarise_data_file(clusters, "cluster") == "3 clusters"
    assert _summarise_data_file(clusters, "cluster_counts") == "3 clusters"

    projection = tmp_path / "pca.dat"
    projection.write_text("\n".join(f"{i} {i * 0.1} {i * 0.2}" for i in range(12)),
                          encoding="utf-8")
    assert _summarise_data_file(projection, "scatter") == "12 frames"

    # Malformed data raises rather than inventing a caption.
    import pytest

    bad = tmp_path / "bad.dat"
    bad.write_text("1\n2\n3\n", encoding="utf-8")
    with pytest.raises(ValueError):
        _summarise_data_file(bad, "scatter")

def test_analysis_view_hosts_sections_and_quick_actions(tmp_path: Path) -> None:
    run = tmp_path / "run"
    run.mkdir()
    server, base_url = start_test_server(run)
    try:
        html = urlopen(f"{base_url}/", timeout=5).read().decode("utf-8")
    finally:
        server.shutdown()
        server.server_close()

    assert 'id="analysis-sections"' in html
    assert 'id="analysis-quick-actions"' in html

    import fastmdxplora.gui as gui_pkg

    script = (Path(gui_pkg.__file__).with_name("static") / "dashboard.js").read_text(
        encoding="utf-8"
    )
    assert "function renderAnalysisSections(" in script
    assert "renderAnalysisSections(payload);" in script


def test_run_title_does_not_repeat_the_product_name(tmp_path: Path) -> None:
    """The sidebar already says FastMDXplora; the run heading names the system."""
    from fastmdxplora.gui.server import _run_title

    assert _run_title(tmp_path, {"system": "1L2Y"}) == "1L2Y"
    assert _run_title(tmp_path, {"system": "/data/protein.pdb"}) == "protein"
    # An explicit report title still wins.
    assert (
        _run_title(tmp_path, {"system": "1L2Y", "options": {"report": {"title": "T4L apo"}}})
        == "T4L apo"
    )


def test_analysis_cards_show_the_publication_figure() -> None:
    """Cards display the figure the analysis wrote, not a second rendering.

    The analyses draw at publication settings; the report's compact variant is
    smaller, so showing it in the card would differ from what opening the
    figure gives you.
    """
    import fastmdxplora.gui as gui_pkg

    script = (Path(gui_pkg.__file__).with_name("static") / "dashboard.js").read_text(
        encoding="utf-8"
    )
    assert "const figure = panel.original_href || panel.href;" in script
    assert "Analysis figure</a>" not in script

    css = (Path(gui_pkg.__file__).with_name("static") / "dashboard.css").read_text(
        encoding="utf-8"
    )
    # A 180px-tall frame scaled a 1950px-wide figure to about 14%, which made
    # axis labels illegible regardless of source DPI.
    assert "aspect-ratio: 6.5 / 4.2;" in css
    assert "minmax(440px, 1fr)" in css


def test_pressure_is_not_advertised_as_a_live_metric(tmp_path: Path) -> None:
    """OpenMM's StateDataReporter cannot report pressure, so no card claims it."""
    run = tmp_path / "run"
    run.mkdir()
    server, base_url = start_test_server(run)
    try:
        html = urlopen(f"{base_url}/", timeout=5).read().decode("utf-8")
    finally:
        server.shutdown()
        server.server_close()
    assert 'data-metric="pressure"' not in html


def test_each_analysis_gets_its_own_section(tmp_path: Path) -> None:
    """Sections follow the analyses, not invented groupings.

    Placement used to depend on figure count: a section holding three figures
    or fewer was merged into "Additional Analysis", so SASA appeared under its
    own heading or a catch-all depending on the run.
    """
    import json

    from fastmdxplora.gui.server import _results_payload

    run = tmp_path / "run"
    (run / "analysis").mkdir(parents=True)
    (run / "manifest.json").write_text(json.dumps({"system": "1L2Y", "phases": []}))
    for folder, figures in {
        "rmsd": ["rmsd"],
        "sasa": ["sasa_total"],
        "cluster": ["kmeans_trajectory", "kmeans_population", "dendrogram"],
    }.items():
        (run / "analysis" / folder).mkdir(parents=True)
        for name in figures:
            (run / "analysis" / folder / f"{name}.png").write_bytes(b"\x89PNG\r\n\x1a\n")
            (run / "analysis" / folder / f"{name}.svg").write_text("<svg/>", encoding="utf-8")

    sections = {s["title"]: s for s in _results_payload(run)["analysis_sections"]}

    # A one-figure analysis keeps its own heading.
    assert "Solvent Accessible Surface Area" in sections
    assert "RMSD" in sections
    assert "Clustering" in sections
    assert "Additional Analysis" not in sections
    assert "Core Metrics" not in sections

    # PNG and SVG of the same plot count once.
    assert len(sections["RMSD"]["panels"]) == 1
    assert len(sections["Clustering"]["panels"]) == 3

    # Anchors are derived from the title, so new analyses need no extra entry.
    assert sections["Solvent Accessible Surface Area"]["anchor"] == (
        "solvent-accessible-surface-area"
    )


def test_long_paths_truncate_instead_of_wrapping() -> None:
    """text-overflow only applies to a single unwrapped line."""
    import fastmdxplora.gui as gui_pkg

    css = (Path(gui_pkg.__file__).with_name("static") / "dashboard.css").read_text(
        encoding="utf-8"
    )
    import re

    for selector in (".footer-value", ".run-id", ".run-title"):
        match = re.search(re.escape(selector) + r"\s*\{([^}]*)\}", css)
        assert match, f"{selector} rule missing"
        body = match.group(1)
        if "text-overflow" in body:
            assert "white-space: nowrap" in body, f"{selector} truncation needs nowrap"

    script = (Path(gui_pkg.__file__).with_name("static") / "dashboard.js").read_text(
        encoding="utf-8"
    )
    # Truncated values stay reachable on hover.
    assert "function setTextWithTooltip(" in script
    assert 'setTextWithTooltip("sidebar-run-name"' in script


def test_summary_card_values_truncate_long_paths() -> None:
    """The Output folder card holds a full path and must not wrap.

    text-overflow needs white-space: nowrap to take effect, and the full value
    is kept in the title attribute so truncation loses nothing.
    """
    import re

    import fastmdxplora.gui as gui_pkg

    css = (Path(gui_pkg.__file__).with_name("static") / "dashboard.css").read_text(
        encoding="utf-8"
    )
    rule = re.search(r"\.metric-card-value\s*\{([^}]*)\}", css)
    assert rule, ".metric-card-value rule missing"
    body = rule.group(1)
    assert "text-overflow: ellipsis" in body
    assert "white-space: nowrap" in body

    script = (Path(gui_pkg.__file__).with_name("static") / "dashboard.js").read_text(
        encoding="utf-8"
    )
    assert 'title="${escapeAttr(value)}"' in script
    assert 'data-kind="${kind}"' in script


class TestTheBrowserCanReachEverySetting:
    """The form was written by hand, one control at a time.

    So it offered eleven of the eighty-three settings that exist, and whole
    phases were missing: there was no way to configure an analysis or a report
    from the browser at all. A field nobody thought to add was simply
    unreachable, and nothing said so.

    The page now reads the schema instead. Adding a field puts a control on the
    page, and there is no second list to keep in step.
    """

    def test_every_schema_field_is_offered(self) -> None:
        from fastmdxplora.config.schema import PHASE_SCHEMAS
        from fastmdxplora.gui.schema_payload import schema_payload

        payload = schema_payload()
        for phase, group in PHASE_SCHEMAS.items():
            offered = {f["name"] for f in payload["phases"][phase]["fields"]}
            declared = set(group.field_names())
            assert offered == declared, (
                f"{phase}: not offered {sorted(declared - offered)}"
            )

    def test_no_phase_is_left_out(self) -> None:
        from fastmdxplora.config.schema import PHASE_SCHEMAS
        from fastmdxplora.gui.schema_payload import schema_payload

        assert set(schema_payload()["phases"]) == set(PHASE_SCHEMAS)

    def test_a_field_carries_what_a_control_needs(self) -> None:
        from fastmdxplora.gui.schema_payload import schema_payload

        payload = schema_payload()
        for phase, block in payload["phases"].items():
            for field in block["fields"]:
                assert field["help"], f"{phase}.{field['name']} has no help"
                assert field["control"], f"{phase}.{field['name']} has no control"

    def test_a_field_with_choices_becomes_a_select(self) -> None:
        from fastmdxplora.gui.schema_payload import schema_payload

        fields = {f["name"]: f for f in
                  schema_payload()["phases"]["setup"]["fields"]}
        assert fields["heterogens"]["control"] == "select"
        assert fields["heterogens"]["choices"] == ["auto", "drop", "keep"]
        assert fields["ph"]["control"] == "number"

    def test_it_is_json(self) -> None:
        """It is sent over the wire, so it must survive the trip."""
        import json

        from fastmdxplora.gui.schema_payload import schema_payload

        restored = json.loads(json.dumps(schema_payload()))
        assert restored["phases"]["analysis"]["fields"]

    def test_the_dashboard_serves_it(self) -> None:
        """Found through the module, and read as the file was written.

        Reading by a path relative to the working directory only works when
        the tests are run from the repository root; and reading without naming
        the encoding uses whatever the machine happens to default to, which is
        ASCII on some and chokes on the first accented character in the file.
        """
        import pathlib

        from fastmdxplora.gui import server

        source = pathlib.Path(server.__file__).read_text(encoding="utf-8")
        assert "/api/schema" in source

    def test_per_analysis_settings_come_from_the_analyses(self) -> None:
        from fastmdxplora.gui.schema_payload import schema_payload

        options = schema_payload()["analysis_options"]
        if not options["available"]:
            import pytest

            pytest.skip(options["reason"])
        assert "cluster" in options["analyses"]
        names = {o["name"] for o in options["analyses"]["cluster"]}
        assert {"n_clusters", "linkage", "features"} <= names

    def test_shared_settings_are_marked_apart(self) -> None:
        """So a form can group what every analysis has separately."""
        from fastmdxplora.gui.schema_payload import schema_payload

        options = schema_payload()["analysis_options"]
        if not options["available"]:
            import pytest

            pytest.skip(options["reason"])
        by_name = {o["name"]: o for o in options["analyses"]["cluster"]}
        assert by_name["n_clusters"]["shared"] is False
        assert by_name["figsize"]["shared"] is True


class TestPointingAtAFolderOfResults:
    """A trajectory that already exists should be enough to start.

    From GROMACS, from AMBER, from someone's own OpenMM script, from a run this
    software did last week. Requiring a simulation to be reproduced here before
    the browser is any use excludes everyone who already has the data, which is
    most of the field.
    """

    @staticmethod
    def _a_finished_run(root):
        (root / "setup").mkdir(parents=True)
        (root / "simulation").mkdir(parents=True)
        (root / "setup" / "system.pdb").write_text("ATOM" * 200)
        (root / "simulation" / "topology.pdb").write_text("ATOM" * 100)
        (root / "simulation" / "production.dcd").write_bytes(b"x" * 5_000_000)
        (root / "simulation" / "equilibration.dcd").write_bytes(b"x" * 100_000)
        (root / "resolved_config.yml").write_text("setup: {}")
        return root

    def test_it_finds_the_trajectory_and_the_topology(self, tmp_path) -> None:
        from fastmdxplora.gui.directory_inspect import inspect_directory

        found = inspect_directory(self._a_finished_run(tmp_path / "run"))
        assert found["ok"] and found["can_analyse"]
        assert {t["name"] for t in found["trajectories"]} == {
            "production.dcd", "equilibration.dcd"}

    def test_it_suggests_the_production_run(self, tmp_path) -> None:
        """The long one, found by size rather than by a naming convention.

        Guessing from the name would only work for runs this software named.
        """
        from fastmdxplora.gui.directory_inspect import inspect_directory

        found = inspect_directory(self._a_finished_run(tmp_path / "run"))
        assert found["suggestion"]["trajectory"].endswith("production.dcd")

    def test_it_pairs_the_topology_beside_it(self, tmp_path) -> None:
        """Not the one in setup/, which describes a different stage."""
        from fastmdxplora.gui.directory_inspect import inspect_directory

        found = inspect_directory(self._a_finished_run(tmp_path / "run"))
        assert found["suggestion"]["topology"].endswith(
            "simulation/topology.pdb")

    def test_it_recognises_a_previous_run(self, tmp_path) -> None:
        from fastmdxplora.gui.directory_inspect import inspect_directory

        found = inspect_directory(self._a_finished_run(tmp_path / "run"))
        assert found["is_previous_run"]
        assert "resolved_config.yml" in found["run_markers"]

    def test_a_lone_structure_is_something_to_set_up(self, tmp_path) -> None:
        """A PDB on its own is a simulation waiting to happen, not a result."""
        from fastmdxplora.gui.directory_inspect import inspect_directory

        (tmp_path / "only.pdb").write_text("ATOM" * 50)
        found = inspect_directory(tmp_path)
        assert found["can_set_up"] and not found["can_analyse"]
        assert found["suggestion"]["trajectory"] is None

    def test_a_missing_folder_says_so(self, tmp_path) -> None:
        from fastmdxplora.gui.directory_inspect import inspect_directory

        found = inspect_directory(tmp_path / "nowhere")
        assert not found["ok"] and "No such directory" in found["error"]

    def test_nothing_opens_the_trajectory(self, tmp_path) -> None:
        """The files here are not trajectories at all -- they are filler.

        Recognising them by extension costs nothing and cannot fail. Opening a
        multi-gigabyte DCD to check would stall the page, and the frame count
        is not needed until an analysis is actually asked for.
        """
        from fastmdxplora.gui.directory_inspect import inspect_directory

        found = inspect_directory(self._a_finished_run(tmp_path / "run"))
        assert found["ok"], "reading the file would have raised here"

    def test_the_dashboard_serves_it(self) -> None:
        import pathlib

        from fastmdxplora.gui import server

        source = pathlib.Path(server.__file__).read_text(encoding="utf-8")
        assert "/api/inspect-directory" in source


class TestTheConfigTheBrowserWritesRunsElsewhere:
    """A laptop is a poor place to run fifty nanoseconds.

    So the browser's job is to decide what to run, and the file it writes has
    to work where the compute is. Which means it must be a file the command
    line accepts -- checked here, with the command line's own validator, while
    the form that produced it is still on the screen, rather than an hour into
    a queue on a cluster.
    """

    @staticmethod
    def _a_reasonable_form():
        return {
            "system": "1ubq",
            "system_id": "ubiquitin",
            "output": "runs/ubq",
            "include": ["setup", "simulation", "analysis"],
            "setup": {"ph": "6.5", "forcefield": "amber14",
                      "keep_water": "true"},
            "simulation": {"temperature_K": "310"},
            "analysis": {"scope": "protein", "include": "rmsd, rmsf"},
        }

    def test_the_command_line_reads_back_what_the_browser_wrote(
        self, tmp_path
    ) -> None:
        from fastmdxplora.config.loader import (
            load_config_file, phase_options, validate_config,
        )
        from fastmdxplora.gui.config_builder import config_yaml

        written = config_yaml(self._a_reasonable_form())
        assert written["ok"], written["error"]

        path = tmp_path / "exploration.yml"
        path.write_text(written["yaml"], encoding="utf-8")

        loaded = load_config_file(path)
        validate_config(loaded)
        options = phase_options(loaded)
        assert options["setup"]["ph"] == 6.5
        assert options["setup"]["forcefield"] == "amber14"
        assert options["analysis"]["scope"] == "protein"

    def test_a_setting_left_alone_is_left_out(self) -> None:
        """A file restating forty defaults buries the three that were chosen.

        Those three are what a reader of the methods section needs.
        """
        from fastmdxplora.gui.config_builder import build_config

        form = self._a_reasonable_form()
        form["setup"]["heterogens"] = "auto"      # the default
        built = build_config(form)
        assert "heterogens" not in built["setup"]
        assert built["setup"]["forcefield"] == "amber14"

    def test_text_from_a_form_becomes_the_right_type(self) -> None:
        """Everything arrives as text, and a type may be several at once."""
        from fastmdxplora.gui.config_builder import build_config

        built = build_config(self._a_reasonable_form())
        assert built["setup"]["ph"] == 6.5
        assert built["setup"]["keep_water"] is True
        assert built["simulation"]["temperature_K"] == 310
        assert built["analysis"]["include"] == ["rmsd", "rmsf"]

    def test_a_setting_that_does_not_exist_is_dropped(self) -> None:
        """The page must not be able to write a file the CLI will refuse."""
        from fastmdxplora.gui.config_builder import build_config

        form = self._a_reasonable_form()
        form["setup"]["not_a_real_setting"] = "x"
        form["not_a_real_phase"] = {"y": 1}
        built = build_config(form)
        assert "not_a_real_setting" not in built["setup"]
        assert "not_a_real_phase" not in built

    def test_a_configuration_that_would_be_refused_says_so_here(self) -> None:
        from fastmdxplora.gui.config_builder import config_yaml

        out = config_yaml({"output": "o", "include": ["not_a_phase"]})
        assert not out["ok"]
        assert out["error"] and out["yaml"] is None

    def test_one_system_is_written_the_way_forty_would_be(self) -> None:
        """So a file for one protein and a file for a batch have one shape."""
        from fastmdxplora.gui.config_builder import build_config

        built = build_config(self._a_reasonable_form())
        assert built["systems"] == [{"system": "1ubq", "id": "ubiquitin"}]

    def test_the_dashboard_serves_it(self) -> None:
        import pathlib

        from fastmdxplora.gui import server

        source = pathlib.Path(server.__file__).read_text(encoding="utf-8")
        assert "/api/config" in source


class TestARunStartsFromTheFileYouCouldTakeAway:
    """One artifact, whether the compute is here or on a cluster.

    A run used to start from a command assembled out of hand-picked flags --
    four setup settings, a few simulation ones -- so the browser could only
    start the kind of run somebody had thought to wire up, and the phase list
    was fixed at setup and simulation before anything else was considered.
    Analysing a trajectory that already existed was not expressible at all.

    Starting from a config settles both: a config names its phases, and the
    file that runs here is the file the download hands over.
    """

    @staticmethod
    def _analyse_what_exists():
        return {
            "output": "runs/analysis",
            "include": ["analysis"],
            "systems": [{"system": "topology.pdb", "id": "existing"}],
            "analysis": {
                "trajectory": "production.dcd",
                "topology": "topology.pdb",
                "scope": "protein",
                "include": "rmsd, rmsf, rg",
            },
        }

    def test_an_analysis_of_existing_data_is_expressible(self, tmp_path) -> None:
        from fastmdxplora.gui.run_from_config import prepare_run

        run = prepare_run(self._analyse_what_exists(), tmp_path)
        assert run["ok"], run["error"]

    def test_what_runs_here_is_what_the_download_gives(self, tmp_path) -> None:
        """Not a second implementation that happens to agree."""
        from fastmdxplora.gui.config_builder import config_yaml
        from fastmdxplora.gui.run_from_config import prepare_run

        state = self._analyse_what_exists()
        run = prepare_run(state, tmp_path)
        assert run["config_yaml"] == config_yaml(state)["yaml"]

    def test_the_config_is_written_beside_the_results(self, tmp_path) -> None:
        """So a run can be repeated from the output directory rather than
        from what somebody remembers typing."""
        import pathlib

        from fastmdxplora.gui.run_from_config import CONFIG_FILENAME, prepare_run

        run = prepare_run(self._analyse_what_exists(), tmp_path)
        written = pathlib.Path(run["config_path"])
        assert written.name == CONFIG_FILENAME
        assert written.parent == tmp_path
        assert written.read_text(encoding="utf-8") == run["config_yaml"]

    def test_the_command_runs_that_file(self, tmp_path) -> None:
        from fastmdxplora.gui.run_from_config import prepare_run

        run = prepare_run(self._analyse_what_exists(), tmp_path)
        assert "--config" in run["command"]
        assert run["config_path"] in run["command"]
        assert "explore" in run["command"]

    def test_the_command_line_accepts_that_file(self, tmp_path) -> None:
        """The parser is asked, rather than the flag being assumed to exist."""
        from fastmdxplora.cli.main import _build_parser
        from fastmdxplora.gui.run_from_config import prepare_run

        run = prepare_run(self._analyse_what_exists(), tmp_path)
        parser = _build_parser()
        # command[3:] drops the interpreter, -m and the module path.
        parsed = parser.parse_args(run["command"][3:])
        assert parsed.config == run["config_path"]

    def test_a_configuration_that_would_be_refused_starts_nothing(
        self, tmp_path
    ) -> None:
        from fastmdxplora.gui.run_from_config import prepare_run

        run = prepare_run({"include": ["not_a_phase"]}, tmp_path)
        assert not run["ok"]
        assert run["command"] is None
        assert not (tmp_path / "exploration.yml").exists()


class TestTheAnalysePageIsWiredToTheSchema:
    """The other way in: data you already have.

    "New Exploration" starts from a protein. This starts from a trajectory --
    from GROMACS, from AMBER, from a run last week -- which is most people, and
    which the dashboard could not do at all.

    Nothing on the page holds a list of what can be set. The controls are drawn
    from what the endpoints report, so an analysis gaining an option puts a
    control on the page. A list somebody has to remember to update is how the
    old form came to offer eleven settings of eighty-three.
    """

    @staticmethod
    def _page():
        import pathlib

        from fastmdxplora.gui import server

        root = pathlib.Path(server.__file__).parent
        return (root / "templates" / "dashboard.html").read_text(encoding="utf-8")

    @staticmethod
    def _script():
        import pathlib

        from fastmdxplora.gui import server

        root = pathlib.Path(server.__file__).parent
        return (root / "static" / "analysis-builder.js").read_text(
            encoding="utf-8")

    def test_the_page_exists_and_is_reachable(self) -> None:
        page = self._page()
        assert 'data-page="analyse"' in page
        assert 'data-view-link="analyse"' in page

    def test_it_does_not_wait_for_a_run(self) -> None:
        """Every other page but the builder is hidden until something runs.

        This one must not be: the whole point is that the run already
        happened, somewhere else.
        """
        import re

        page = self._page()
        link = re.search(r'<a[^>]*data-view-link="analyse"[^>]*>', page)
        assert link, "no navigation link for the analyse page"
        assert "data-requires-run" not in link.group(0)

    def test_the_script_is_loaded(self) -> None:
        assert "/static/analysis-builder.js" in self._page()

    def test_it_asks_the_endpoints_rather_than_holding_a_list(self) -> None:
        script = self._script()
        for endpoint in ("/api/schema", "/api/inspect-directory", "/api/config"):
            assert endpoint in script, f"{endpoint} is never called"

    def test_no_analysis_is_named_in_the_page_or_the_script(self) -> None:
        """The moment one is, the page has a list to keep in step.

        Every analysis on the page comes from /api/schema, which reads them
        from the analyses themselves.
        """
        import fastmdxplora.analysis  # noqa: F401
        from fastmdxplora.analysis.orchestrator import available_analyses

        script = self._script()
        hardcoded = [
            name for name in available_analyses()
            if f'"{name}"' in script or f"'{name}'" in script
        ]
        assert not hardcoded, (
            f"these analyses are written into the page: {hardcoded}"
        )

    def test_every_control_it_draws_has_somewhere_to_get_its_values(
        self,
    ) -> None:
        """Each control the script knows how to draw must be a control the
        payload actually asks for, and the other way about."""
        from fastmdxplora.gui.schema_payload import schema_payload

        script = self._script()
        payload = schema_payload()
        wanted = {
            field["control"]
            for block in payload["phases"].values()
            for field in block["fields"]
        }
        # These are the ones this page draws; the rest belong to the setup and
        # simulation forms, which are still written by hand.
        drawable = {"select", "number", "text", "checkbox"}
        assert drawable <= wanted, "the page draws a control nothing asks for"
        for control in drawable:
            assert control in script, f"the script cannot draw a {control}"

    def test_the_page_offers_the_scopes_the_schema_declares(self) -> None:
        from fastmdxplora.config.schema import PHASE_SCHEMAS

        page = self._page()
        scope = next(f for f in PHASE_SCHEMAS["analysis"].fields
                     if f.name == "scope")
        for choice in scope.choices:
            assert f'value="{choice}"' in page, f"scope {choice} is not offered"

    def test_the_styles_it_uses_are_defined(self) -> None:
        import pathlib

        from fastmdxplora.gui import server

        css = (pathlib.Path(server.__file__).parent / "static"
               / "dashboard.css").read_text(encoding="utf-8")
        for name in ("analyse-choices", "analyse-row", "analyse-toggle",
                     "analyse-settings"):
            assert f".{name}" in css, f"{name} is used but never styled"

    def test_reduced_motion_is_respected(self) -> None:
        import pathlib

        from fastmdxplora.gui import server

        css = (pathlib.Path(server.__file__).parent / "static"
               / "dashboard.css").read_text(encoding="utf-8")
        assert "prefers-reduced-motion" in css


class TestEveryPageCanBeReached:
    """A section in the document is not the same as a page you can get to.

    The router kept its own list of page names and sent anything not on it to
    the overview. So adding a section was not enough: the new page's link
    appeared, was clickable, and quietly showed the overview instead -- which
    looks exactly like a link wired to the wrong place.

    The pages are now read from the document, so a section that exists is a
    section you can reach.
    """

    @staticmethod
    def _files():
        import pathlib

        from fastmdxplora.gui import server

        root = pathlib.Path(server.__file__).parent
        return (
            (root / "templates" / "dashboard.html").read_text(encoding="utf-8"),
            (root / "static" / "dashboard.js").read_text(encoding="utf-8"),
        )

    def test_the_router_does_not_keep_its_own_list(self) -> None:
        import re

        _page, script = self._files()
        declared = re.search(r"pages:\s*\[([^\]]*)\]", script)
        assert declared, "the page list is gone entirely; check this test"
        assert not declared.group(1).strip(), (
            "the router names its pages, so a new section will not be reachable "
            f"until somebody remembers to add it here: {declared.group(1)}"
        )

    def test_the_router_reads_them_from_the_document(self) -> None:
        _page, script = self._files()
        assert 'state.pages = $$(\'.page\')' in script

    def test_every_section_has_a_link_and_every_link_a_section(self) -> None:
        import re

        page, _script = self._files()
        sections = set(re.findall(r'<section class="page" data-page="([^"]+)"', page))
        links = set(re.findall(r'data-view-link="([^"]+)"', page))
        # Links may point at a page from elsewhere on the page, but a section
        # nobody links to cannot be opened and a link to nothing is a dead end.
        assert sections, "no pages found at all"
        assert links <= sections, f"links with no section: {sorted(links - sections)}"
        assert sections <= links, f"sections with no link: {sorted(sections - links)}"

    def test_the_analyse_page_is_among_them(self) -> None:
        page, _script = self._files()
        assert 'data-page="analyse"' in page
        assert 'data-view-link="analyse"' in page


class TestChoosingAFolderWithoutTypingOne:
    """A browser cannot offer a dialog for a folder the server will read.

    The one it has uploads files to a page; this needs a path the process can
    open. So the walking happens here and the page draws it.
    """

    def test_it_lists_folders_and_the_way_back(self, tmp_path) -> None:
        from fastmdxplora.gui.browse import browse

        (tmp_path / "runs").mkdir()
        (tmp_path / "notes").mkdir()
        listing = browse(tmp_path)
        assert listing["ok"]
        assert {e["name"] for e in listing["entries"]} == {"runs", "notes"}
        assert listing["parent"] == str(tmp_path.parent)

    def test_it_says_which_folders_hold_data(self, tmp_path) -> None:
        """So a list of run folders can be read without opening each one."""
        from fastmdxplora.gui.browse import browse

        (tmp_path / "has_data").mkdir()
        (tmp_path / "has_data" / "run.dcd").write_bytes(b"x")
        (tmp_path / "has_data" / "top.pdb").write_text("ATOM")
        (tmp_path / "empty").mkdir()

        by_name = {e["name"]: e for e in browse(tmp_path)["entries"]}
        assert by_name["has_data"]["trajectories"] == 1
        assert by_name["has_data"]["structures"] == 1
        assert by_name["empty"]["trajectories"] == 0

    def test_nothing_is_opened_to_count_them(self, tmp_path) -> None:
        """The files here are not trajectories. Counting by extension cannot
        fail; reading them would."""
        from fastmdxplora.gui.browse import browse

        (tmp_path / "run").mkdir()
        (tmp_path / "run" / "not_really.dcd").write_text("this is not a DCD")
        assert browse(tmp_path)["entries"][0]["trajectories"] == 1

    def test_hidden_folders_are_left_out(self, tmp_path) -> None:
        from fastmdxplora.gui.browse import browse

        (tmp_path / ".git").mkdir()
        (tmp_path / "runs").mkdir()
        assert [e["name"] for e in browse(tmp_path)["entries"]] == ["runs"]

    def test_a_file_means_the_folder_holding_it(self, tmp_path) -> None:
        """Pasting a path from a terminal often lands on a file."""
        from fastmdxplora.gui.browse import browse

        target = tmp_path / "topology.pdb"
        target.write_text("ATOM")
        assert browse(target)["path"] == str(tmp_path.resolve())

    def test_no_path_starts_somewhere_readable(self) -> None:
        import pathlib

        from fastmdxplora.gui.browse import browse

        assert browse()["path"] == str(pathlib.Path.home().resolve())

    def test_a_folder_that_is_not_there_says_so(self, tmp_path) -> None:
        from fastmdxplora.gui.browse import browse

        listing = browse(tmp_path / "nowhere")
        assert not listing["ok"] and "No such folder" in listing["error"]

    def test_the_dashboard_serves_it(self) -> None:
        import pathlib

        from fastmdxplora.gui import server

        source = pathlib.Path(server.__file__).read_text(encoding="utf-8")
        assert "/api/browse" in source

    def test_both_folder_fields_can_be_browsed(self) -> None:
        import pathlib

        from fastmdxplora.gui import server

        page = (pathlib.Path(server.__file__).parent / "templates"
                / "dashboard.html").read_text(encoding="utf-8")
        assert 'id="analyse-browse"' in page
        assert 'id="analyse-browse-output"' in page


class TestASettingShowsWhatItWillDo:
    """A control that cannot say what happens if left alone is not much help.

    Two settings hid their default behind None and filled it in afterwards, so
    anything reading the signature -- a form, a config template, the help --
    saw no default at all. Clustering was worse than silent: its docstring
    named one method where the code ran two.
    """

    def test_clustering_declares_the_methods_it_runs(self) -> None:
        from fastmdxplora.analysis.cluster import Cluster

        assert Cluster().methods == ["kmeans", "hierarchical"]

    def test_and_the_description_agrees(self) -> None:
        import fastmdxplora.analysis  # noqa: F401
        from fastmdxplora.analysis.describe import describe_analysis

        option = next(o for o in describe_analysis("cluster")
                      if o.name == "methods")
        assert option.default == ("kmeans", "hierarchical")

    def test_dimensionality_reduction_too(self) -> None:
        import fastmdxplora.analysis  # noqa: F401
        from fastmdxplora.analysis.describe import describe_analysis
        from fastmdxplora.analysis.dimred import DimRed

        option = next(o for o in describe_analysis("dimred")
                      if o.name == "methods")
        assert option.default == ("pca",)
        assert DimRed().methods == ["pca"]

    def test_only_the_genuinely_undecidable_have_no_default(self) -> None:
        """A ligand's residue name depends on the structure, so it has none
        until one is read. Everything else should be able to say."""
        import fastmdxplora.analysis  # noqa: F401
        from fastmdxplora.analysis.describe import describe_all

        silent = {
            f"{name}.{o.name}"
            for name, options in describe_all().items()
            for o in options
            if o.default is None and o.owner != "Analysis"
            and o.name != "ligand_resname"
        }
        assert not silent, f"these cannot say what they would do: {sorted(silent)}"


class TestAnAnalysisExplainsItselfOnThePage:
    """Every analysis module opens by saying what its measurement means.

    What a rising profile indicates, which paper the definition comes from,
    why the default is the default. That was written for somebody reading the
    source, and it is exactly what somebody choosing between checkboxes needs.
    It was not on the page.
    """

    def test_every_analysis_has_something_to_say(self) -> None:
        import fastmdxplora.analysis  # noqa: F401
        from fastmdxplora.analysis.describe import explain_analysis
        from fastmdxplora.analysis.orchestrator import available_analyses

        silent = []
        for name in available_analyses():
            explanation = explain_analysis(name)
            if not explanation["summary"] or not explanation["detail"]:
                silent.append(name)
        assert not silent, f"these would appear with no explanation: {silent}"

    def test_the_explanation_reaches_the_payload(self) -> None:
        from fastmdxplora.gui.schema_payload import schema_payload

        options = schema_payload()["analysis_options"]
        if not options["available"]:
            import pytest

            pytest.skip(options["reason"])
        qvalue = options["explanations"]["qvalue"]
        assert "folding-state" in qvalue["detail"], qvalue["detail"][:80]
        assert qvalue["title"] == "Native contact fraction (Q)"

    def test_the_page_draws_it(self) -> None:
        import pathlib

        from fastmdxplora.gui import server

        root = pathlib.Path(server.__file__).parent
        script = (root / "static" / "analysis-builder.js").read_text(
            encoding="utf-8")
        css = (root / "static" / "dashboard.css").read_text(encoding="utf-8")
        assert "explanations" in script
        assert "analyse-explains" in script
        assert ".analyse-explains" in css


class TestChoosingBetweenSeveralTrajectories:
    """A run with replicates has a production.dcd in each of them.

    Listing them by name alone is no choice at all.
    """

    @staticmethod
    def _replicates(root):
        # Sizes must actually differ for "the longest is chosen" to mean
        # anything; naming them rep1..rep3 makes len() the same for each.
        for index, name in enumerate(("rep1", "rep2", "rep3"), start=1):
            (root / name).mkdir(parents=True)
            (root / name / "production.dcd").write_bytes(b"x" * (1000 * index))
            (root / name / "topology.pdb").write_text("ATOM")
        return root

    def test_all_of_them_are_offered(self, tmp_path) -> None:
        from fastmdxplora.gui.directory_inspect import inspect_directory

        found = inspect_directory(self._replicates(tmp_path))
        assert len(found["trajectories"]) == 3

    def test_they_are_told_apart_by_where_they_sit(self, tmp_path) -> None:
        """The paths differ even where the names do not, and the page labels
        each by its path below the folder that was searched."""
        import pathlib

        from fastmdxplora.gui import server
        from fastmdxplora.gui.directory_inspect import inspect_directory

        found = inspect_directory(self._replicates(tmp_path))
        names = {pathlib.Path(t["path"]).name for t in found["trajectories"]}
        paths = {t["path"] for t in found["trajectories"]}
        assert names == {"production.dcd"}, "the fixture should share a name"
        assert len(paths) == 3, "and differ by path"

        script = (pathlib.Path(server.__file__).parent / "static"
                  / "analysis-builder.js").read_text(encoding="utf-8")
        assert "entry.path" in script, "the page labels by name alone"

    def test_the_largest_is_chosen_but_not_forced(self, tmp_path) -> None:
        from fastmdxplora.gui.directory_inspect import inspect_directory

        found = inspect_directory(self._replicates(tmp_path))
        assert found["suggestion"]["trajectory"].endswith("rep3/production.dcd")
        # and the others remain on offer
        assert len(found["trajectories"]) == 3


class TestFoldersBelowTheOneYouAreLookingAt:
    """A run keeps its trajectory in simulation/, a package keeps its data in
    <name>/data/, and the folder somebody is looking at is the one above.

    Counting only what sits directly inside marked almost nothing.
    """

    def test_data_a_few_levels_down_is_found(self, tmp_path) -> None:
        from fastmdxplora.gui.browse import browse

        deep = tmp_path / "_v1" / "fastmdanalysis" / "data"
        deep.mkdir(parents=True)
        (deep / "trp_cage.dcd").write_bytes(b"x")
        (deep / "trp_cage.pdb").write_text("ATOM")

        entry = next(e for e in browse(tmp_path)["entries"] if e["name"] == "_v1")
        assert entry["trajectories"] == 1
        assert entry["structures"] == 1

    def test_it_does_not_walk_the_whole_disk(self, tmp_path) -> None:
        """Bounded, because this runs once per row of a listing."""
        from fastmdxplora.gui.browse import browse

        far = tmp_path / "archive" / "a" / "b" / "c" / "d"
        far.mkdir(parents=True)
        (far / "buried.dcd").write_bytes(b"x")

        entry = next(e for e in browse(tmp_path)["entries"]
                     if e["name"] == "archive")
        assert entry["trajectories"] == 0

    def test_it_stops_counting_once_the_answer_is_clear(self, tmp_path) -> None:
        from fastmdxplora.gui.browse import browse

        many = tmp_path / "many"
        many.mkdir()
        for i in range(60):
            (many / f"frame{i}.dcd").write_bytes(b"x")

        entry = next(e for e in browse(tmp_path)["entries"] if e["name"] == "many")
        assert entry["trajectories"] < 60
        assert entry["truncated"]


class TestBothKindsOfConfigFile:
    """Two files answer two questions, and both run.

    Left to itself the file names only what was decided, because a file
    restating forty settings the software would have chosen anyway hides the
    three that were chosen. Asked for in full it names everything, which is
    the file to keep beside a result: a default that moves in a later version
    cannot then change what the file meant.
    """

    @staticmethod
    def _state():
        return {
            "system": "1ubq",
            "output": "runs/x",
            "include": ["setup", "simulation", "analysis"],
            "analysis": {"scope": "protein"},
        }

    def test_the_short_one_names_only_what_changed(self) -> None:
        from fastmdxplora.gui.config_builder import build_config

        built = build_config(self._state())
        assert built["analysis"] == {"scope": "protein"}
        assert "setup" not in built

    def test_the_full_one_names_every_phase_that_runs(self) -> None:
        """A phase the form never touched still runs, and still uses values
        worth recording."""
        from fastmdxplora.gui.config_builder import build_config

        built = build_config(self._state(), full=True)
        for phase in ("setup", "simulation", "analysis"):
            assert built[phase], f"{phase} runs but is not recorded"
        assert built["setup"]["ph"] == 7.4
        assert built["analysis"]["scope"] == "protein"

    def test_a_phase_that_does_not_run_is_not_recorded(self) -> None:
        from fastmdxplora.gui.config_builder import build_config

        built = build_config(self._state(), full=True)
        assert "report" not in built

    def test_both_are_accepted_by_the_command_line(self, tmp_path) -> None:
        from fastmdxplora.config.loader import load_config_file, validate_config
        from fastmdxplora.gui.config_builder import config_yaml

        for index, full in enumerate((False, True)):
            written = config_yaml(self._state(), full=full)
            assert written["ok"], written["error"]
            path = tmp_path / f"config{index}.yml"
            path.write_text(written["yaml"], encoding="utf-8")
            validate_config(load_config_file(path))

    def test_they_agree_on_what_the_run_will_do(self, tmp_path) -> None:
        """The full file must not change the run, only describe it.

        Every value it adds is the value the software would have used anyway,
        so the two files have to resolve to the same settings.
        """
        from fastmdxplora.config.schema import PHASE_SCHEMAS
        from fastmdxplora.gui.config_builder import build_config

        short = build_config(self._state())
        full = build_config(self._state(), full=True)
        for phase, group in PHASE_SCHEMAS.items():
            defaults = {f.name: f.default for f in group.fields}
            for name, value in full.get(phase, {}).items():
                if name == "options":
                    # The nested per-analysis settings have no schema default
                    # to compare against; they are checked against the
                    # analyses themselves in TestTheFullConfigReachesTheNested-
                    # Settings, which is where they come from.
                    continue
                chosen = short.get(phase, {}).get(name)
                assert value == (chosen if chosen is not None else defaults[name]), (
                    f"{phase}.{name} differs between the two files"
                )

    def test_the_header_says_which_kind_it_is(self) -> None:
        from fastmdxplora.gui.config_builder import config_yaml

        assert "Only settings that differ" in config_yaml(self._state())["yaml"]
        assert "Every setting is named" in config_yaml(
            self._state(), full=True)["yaml"]

    def test_the_page_offers_the_choice(self) -> None:
        import pathlib

        from fastmdxplora.gui import server

        root = pathlib.Path(server.__file__).parent
        page = (root / "templates" / "dashboard.html").read_text(encoding="utf-8")
        script = (root / "static" / "analysis-builder.js").read_text(
            encoding="utf-8")
        assert 'id="analyse-full-config"' in page
        assert "body.full" in script


class TestItIsCalledTheGUI:
    """A dashboard is the live view inside it, not the thing itself."""

    def test_the_page_is_titled_that_way(self) -> None:
        import pathlib

        from fastmdxplora.gui import server

        page = (pathlib.Path(server.__file__).parent / "templates"
                / "dashboard.html").read_text(encoding="utf-8")
        assert "<title>FastMDXplora GUI</title>" in page

    def test_the_config_it_writes_says_so(self) -> None:
        from fastmdxplora.gui.config_builder import config_yaml

        written = config_yaml({"output": "o", "systems": [{"system": "1ubq"}]})
        assert "FastMDXplora GUI" in written["yaml"]
        assert "dashboard" not in written["yaml"].lower()


class TestTheFullConfigReachesTheNestedSettings:
    """"Every setting" has to mean the nested ones too.

    A full config that stopped at the phase boundary recorded the scope and
    the stride and said nothing about how the clustering was done -- which is
    where the decisions that change a result actually live.
    """

    @staticmethod
    def _clustering():
        return {
            "system": "x", "output": "o", "include": ["analysis"],
            "analysis": {"scope": "protein", "include": "cluster"},
        }

    def test_the_analysis_settings_are_recorded(self) -> None:
        from fastmdxplora.gui.config_builder import build_config

        built = build_config(self._clustering(), full=True)
        cluster = built["analysis"]["options"]["cluster"]
        assert {"n_clusters", "linkage", "features", "methods"} <= set(cluster)

    def test_they_are_the_values_the_run_would_use(self) -> None:
        import fastmdxplora.analysis  # noqa: F401
        from fastmdxplora.analysis.cluster import Cluster
        from fastmdxplora.gui.config_builder import build_config

        recorded = build_config(self._clustering(), full=True)
        cluster = recorded["analysis"]["options"]["cluster"]
        running = Cluster()
        assert cluster["n_clusters"] == running.n_clusters
        assert cluster["linkage"] == running.linkage
        assert cluster["features"] == running.features
        assert cluster["methods"] == running.methods

    def test_a_choice_still_wins_over_the_default(self) -> None:
        from fastmdxplora.gui.config_builder import build_config

        state = self._clustering()
        state["analysis"]["options"] = {"cluster": {"n_clusters": "8"}}
        cluster = build_config(state, full=True)["analysis"]["options"]["cluster"]
        assert cluster["n_clusters"] == 8
        assert cluster["linkage"] == "average"     # untouched, still recorded

    def test_a_number_from_a_form_is_recorded_as_a_number(self) -> None:
        """Left as text the file would say n_clusters: '8', which runs -- the
        analysis casts it -- but reads as though somebody meant the string."""
        from fastmdxplora.gui.config_builder import build_config

        state = self._clustering()
        state["analysis"]["options"] = {"cluster": {"n_clusters": "8",
                                                    "eps": "0.35"}}
        cluster = build_config(state, full=True)["analysis"]["options"]["cluster"]
        assert isinstance(cluster["n_clusters"], int)
        assert isinstance(cluster["eps"], float)

    def test_only_the_chosen_analyses_are_recorded(self) -> None:
        from fastmdxplora.gui.config_builder import build_config

        built = build_config(self._clustering(), full=True)
        assert set(built["analysis"]["options"]) == {"cluster"}

    def test_the_short_config_leaves_them_out(self) -> None:
        from fastmdxplora.gui.config_builder import build_config

        built = build_config(self._clustering())
        assert "options" not in built["analysis"]

    def test_it_still_validates(self, tmp_path) -> None:
        from fastmdxplora.config.loader import load_config_file, validate_config
        from fastmdxplora.gui.config_builder import config_yaml

        written = config_yaml(self._clustering(), full=True)
        assert written["ok"], written["error"]
        path = tmp_path / "full.yml"
        path.write_text(written["yaml"], encoding="utf-8")
        validate_config(load_config_file(path))


class TestWhatTheGUIOpensOn:
    """With nothing running, the thing to do is start something.

    The overview was the only page left visible in the markup, so it was what
    the GUI opened on -- an overview of a run, before any state had loaded and
    usually when there was no run at all.
    """

    @staticmethod
    def _page():
        import pathlib

        from fastmdxplora.gui import server

        return (pathlib.Path(server.__file__).parent / "templates"
                / "dashboard.html").read_text(encoding="utf-8")

    def test_only_one_page_is_visible_before_anything_loads(self) -> None:
        import re

        sections = re.findall(
            r'<section class="page" data-page="([^"]+)"([^>]*)>', self._page())
        visible = [name for name, attrs in sections if "hidden" not in attrs]
        assert visible == ["builder"], f"visible at load: {visible}"

    def test_the_router_decides_once_it_knows(self) -> None:
        import pathlib

        from fastmdxplora.gui import server

        script = (pathlib.Path(server.__file__).parent / "static"
                  / "dashboard.js").read_text(encoding="utf-8")
        # A run to look at, a page asked for by name, or nothing yet.
        assert 'else if (activeRun) navigate("overview");' in script
        assert 'else navigate("builder");' in script


class TestOpeningTheGUIWithNothingRunning:
    """`fastmdx gui` with no --output has no run to watch.

    The working directory is merely where the command was typed. Passing it as
    the active run made the GUI report one, and the router quite correctly
    showed the overview of it -- an overview of nothing. Marking the overview
    hidden did not help, because the router put it back a moment later.
    """

    def test_the_command_says_there_is_no_run(self, tmp_path, monkeypatch) -> None:
        import fastmdxplora.gui.server as server_module
        from fastmdxplora.cli.main import _build_parser

        seen = {}

        def _capture(**kwargs):
            seen.update(kwargs)

        monkeypatch.setattr(server_module, "serve_dashboard", _capture)
        monkeypatch.chdir(tmp_path)

        # `fastmdxplora.cli.main` is a function re-exported by the package,
        # so importing it as a module gets the function instead.
        from fastmdxplora.cli.main import _cmd_gui

        args = _build_parser().parse_args(["gui", "--no-browser"])
        _cmd_gui(args)
        assert seen["home_mode"] is True, (
            "the working directory was passed as an active run"
        )

    def test_an_output_given_is_a_run_to_watch(self, tmp_path, monkeypatch) -> None:
        import fastmdxplora.gui.server as server_module
        from fastmdxplora.cli.main import _build_parser

        seen = {}
        monkeypatch.setattr(server_module, "serve_dashboard",
                            lambda **kw: seen.update(kw))

        from fastmdxplora.cli.main import _cmd_gui

        args = _build_parser().parse_args(
            ["gui", "--no-browser", "--output", str(tmp_path)])
        _cmd_gui(args)
        assert seen["home_mode"] is False
        assert seen["output"] == tmp_path

    def test_the_runtime_reports_no_run_in_that_mode(self, tmp_path) -> None:
        from fastmdxplora.gui.server import DashboardRuntime

        runtime = DashboardRuntime(
            workspace_root=tmp_path, exploration_root=tmp_path.parent,
            active_root=None,
        )
        assert runtime.snapshot()["active_run"] is None

    def test_and_reports_one_when_there_is_one(self, tmp_path) -> None:
        from fastmdxplora.gui.server import DashboardRuntime

        runtime = DashboardRuntime(
            workspace_root=tmp_path, exploration_root=tmp_path.parent,
            active_root=tmp_path,
        )
        assert runtime.snapshot()["active_run"] == str(tmp_path.resolve())
