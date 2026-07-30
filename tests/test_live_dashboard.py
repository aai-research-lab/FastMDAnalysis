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
from fastmdxplora.live import protein_preview
from fastmdxplora.live.ligand_detection import detect_ligands, normalise_ligand_resname
from fastmdxplora.live.live_frames import (
    dashboard_display_pdb,
    read_live_frame_history,
    write_live_frame,
    write_openmm_live_frame,
)
from fastmdxplora.live.protein_preview import protein_preview_payload
from fastmdxplora.live.server import (
    DashboardConfig,
    make_handler,
    start_dashboard_session,
    start_test_server,
)
from fastmdxplora.live.structure_info import count_structure, ligand_atom_counts
from fastmdxplora.simulation.runner import _maybe_write_live_frame
from fastmdxplora.live.telemetry import (
    TelemetryWriter,
    analyze_health,
    read_events,
    read_metrics,
    read_status,
)
from fastmdxplora.live.trajectory_playback import (
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


def test_dashboard_serve_command_parses() -> None:
    args = _build_parser().parse_args(
        ["dashboard", "serve", "--output", "local_runs/my_run", "--port", "8766"]
    )
    assert args.command == "dashboard"
    assert args.dashboard_command == "serve"
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
    import fastmdxplora.live.trajectory_playback as mod
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
    import fastmdxplora.live as live_pkg
    layout = Path(live_pkg.__file__).with_name("templates") / "dashboard.html"
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
    assert "AAI Research Lab" in html or "aai-research-logo" in html
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
    assert "/static/aai-research-logo.svg" in html or "/static/aai-research-mark.svg" in html

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
    css_path = Path(make_handler.__module__.split(".")[0])
    import fastmdxplora.live as live_pkg
    css = (Path(live_pkg.__file__).with_name("static") / "dashboard.css").read_text(encoding="utf-8")
    for token in ("--background-primary", "--accent-cyan", "--accent-violet",
                  "--text-primary", "prefers-reduced-motion", "Inter"):
        assert token in css


def test_dashboard_endpoint_no_inline_html_css() -> None:
    """The server.py module no longer contains the giant HTML f-string.

    We assert the absence of known-fingerprinted fragments from the old
    inline template — a tautology-free check.
    """
    server_path = Path(make_handler.__module__.split(".")[0])
    import fastmdxplora.live as live_pkg
    src = (Path(live_pkg.__file__).with_name("server.py")).read_text(encoding="utf-8")
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


def test_aai_logo_assets_are_packaged(tmp_path: Path) -> None:
    import fastmdxplora.live as live_pkg
    static = Path(live_pkg.__file__).with_name("static")
    assert (static / "aai-research-logo.svg").is_file()
    assert (static / "aai-research-mark.svg").is_file()
    assert (static / "dashboard.css").is_file()
    assert (static / "dashboard.js").is_file()
    assert (static / "charts.js").is_file()
    assert (static / "molecule-viewer.js").is_file()


def test_static_logo_endpoint_serves_svg(tmp_path: Path) -> None:
    run = tmp_path / "run"
    run.mkdir()
    server, base_url = start_test_server(run)
    try:
        response = urlopen(f"{base_url}/static/aai-research-logo.svg", timeout=5)
        body = response.read().decode("utf-8", errors="ignore")
    finally:
        server.shutdown()
        server.server_close()
    assert response.headers["Content-Type"].startswith("image/svg+xml") or "svg" in body.lower()
    assert "AAI" in body

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
    import fastmdxplora.live as live_pkg

    root = Path(live_pkg.__file__).parent
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

    import fastmdxplora.live.live_frames as live_frames_module

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
    html = (root / "src" / "fastmdxplora" / "live" / "templates" / "dashboard.html").read_text(encoding="utf-8")
    dashboard_js = (root / "src" / "fastmdxplora" / "live" / "static" / "dashboard.js").read_text(encoding="utf-8")
    viewer_js = (root / "src" / "fastmdxplora" / "live" / "static" / "molecule-viewer.js").read_text(encoding="utf-8")

    assert 'id="download-all-svg"' in html
    assert "Download SVG" in dashboard_js
    assert "const hadModel" in viewer_js
    assert "opts.center !== false || !hadModel" in viewer_js
    assert "resolveProteinSelection" in viewer_js
