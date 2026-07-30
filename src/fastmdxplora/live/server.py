"""Dependency-free localhost dashboard server.

The HTML shell, CSS, and JavaScript live in
:mod:`fastmdxplora.live.templates` and :mod:`fastmdxplora.live.static`.
``_dashboard_shell`` reads ``dashboard.html`` from disk on first request
and caches the result.
"""

from __future__ import annotations

import io
import json
import logging
import mimetypes
import os
import subprocess
import sys
import threading
import time
import zipfile
from dataclasses import dataclass
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Any
from urllib.parse import unquote, urlparse

from fastmdxplora.live.launcher import (
    DashboardRuntime,
    build_launcher_command,
    launcher_defaults,
    validate_launcher_payload,
)
from fastmdxplora.live.ligand_detection import detect_ligands, normalise_ligand_resname
from fastmdxplora.live.live_frames import live_frame_exists, read_live_frame_index
from fastmdxplora.live.protein_preview import find_structure, protein_preview_payload
from fastmdxplora.live.structure_info import count_structure, ligand_atom_counts
from fastmdxplora.live.telemetry import analyze_health, read_events, read_metrics, read_status
from fastmdxplora.live.trajectory_playback import playback_info

logger = logging.getLogger("fastmdxplora.live.server")

PLOT_TITLE_ALIASES = {
    "rmsd": "RMSD",
    "rmsf": "RMSF",
    "rg": "Radius of gyration",
    "hbonds": "Hydrogen bonds",
    "hbond": "Hydrogen bonds",
    "sasa": "SASA",
    "ss": "Secondary structure",
    "secondary_structure": "Secondary structure",
    "dimred_pca": "PCA",
    "pca": "PCA",
    "dimred_mds": "MDS",
    "mds": "MDS",
    "dimred_tsne": "t-SNE",
    "tsne": "t-SNE",
    "cluster_kmeans": "KMeans trajectory scatter",
    "kmeans_trajectory": "KMeans trajectory scatter",
    "kmeans_population": "KMeans population",
    "cluster_kmeans_counts": "KMeans population",
    "cluster_hierarchical": "Hierarchical trajectory scatter",
    "hierarchical_trajectory": "Hierarchical trajectory scatter",
    "hierarchical_population": "Hierarchical population",
    "cluster_hierarchical_counts": "Hierarchical population",
    "cluster_hierarchical_dendrogram": "Hierarchical dendrogram",
    "hierarchical_dendrogram": "Hierarchical dendrogram",
    "cluster_dbscan": "DBSCAN trajectory scatter",
    "dbscan_trajectory": "DBSCAN trajectory scatter",
    "dbscan_population": "DBSCAN population",
    "cluster_dbscan_counts": "DBSCAN population",
    "qvalue": "Native contacts",
    "native_contacts": "Native contacts",
    "dihedrals": "Dihedrals",
    "ramachandran": "Ramachandran plot",
}

DASHBOARD_ASSET_TITLES = {
    "rmsd_dashboard": "RMSD",
    "rmsf_dashboard": "RMSF",
    "rg_dashboard": "Radius of gyration",
    "hbonds_dashboard": "Hydrogen bonds",
    "sasa_dashboard": "SASA",
    "pca_dashboard": "PCA",
    "mds_dashboard": "MDS",
    "tsne_dashboard": "t-SNE",
    "kmeans_trajectory_dashboard": "KMeans trajectory scatter",
    "kmeans_population_dashboard": "KMeans population",
    "hierarchical_trajectory_dashboard": "Hierarchical trajectory scatter",
    "hierarchical_population_dashboard": "Hierarchical population",
    "hierarchical_dendrogram_dashboard": "Hierarchical dendrogram",
    "cluster_hierarchical_dendrogram_dashboard": "Hierarchical dendrogram",
    "dbscan_trajectory_dashboard": "DBSCAN trajectory scatter",
    "dbscan_population_dashboard": "DBSCAN population",
    "ss_dashboard": "Secondary structure",
    "qvalue_dashboard": "Native contacts",
    "dihedrals_dashboard": "Dihedrals",
}

PLOT_CATEGORY_BY_TITLE = {
    "RMSD": "Core Metrics",
    "RMSF": "Core Metrics",
    "Radius of gyration": "Core Metrics",
    "Hydrogen bonds": "Core Metrics",
    "KMeans trajectory scatter": "Clustering",
    "KMeans population": "Clustering",
    "Hierarchical trajectory scatter": "Clustering",
    "Hierarchical population": "Clustering",
    "Hierarchical dendrogram": "Clustering",
    "DBSCAN trajectory scatter": "Clustering",
    "DBSCAN population": "Clustering",
    "PCA": "Additional Analysis",
    "MDS": "Additional Analysis",
    "t-SNE": "Additional Analysis",
    "SASA": "Additional Analysis",
    "Secondary structure": "Additional Analysis",
    "Native contacts": "Additional Analysis",
    "Dihedrals": "Additional Analysis",
    "Ramachandran plot": "Additional Analysis",
}

KEY_PLOT_TITLES = {"RMSD", "RMSF", "Radius of gyration", "Hydrogen bonds", "PCA", "SASA"}

DASHBOARD_TEMPLATE_PATH = Path(__file__).with_name("templates") / "dashboard.html"


@dataclass
class DashboardConfig:
    """Per-dashboard-server knobs supplied by the CLI."""

    ligand_resname: str | None = None
    include_cofactors: bool = False
    binding_pocket_cutoff_A: float = 5.0
    max_browser_frames: int = 200
    refresh_seconds: float = 3.0

    @property
    def binding_pocket_cutoff_m(self) -> float:
        return float(self.binding_pocket_cutoff_A)


@dataclass
class DashboardSession:
    """Background local dashboard server session."""

    server: ThreadingHTTPServer
    thread: threading.Thread
    root: Path
    host: str
    port: int
    requested_port: int
    config: DashboardConfig | None = None
    runtime: DashboardRuntime | None = None

    @property
    def url(self) -> str:
        return f"http://{self.host}:{self.port}"

    @property
    def port_was_changed(self) -> bool:
        # Port 0 explicitly asks the operating system for any free port; that
        # is not a fallback caused by a conflict and should not be reported as
        # one to the user.
        return self.requested_port != 0 and self.port != self.requested_port

    def stop(self) -> None:
        self.server.shutdown()
        self.server.server_close()
        self.thread.join(timeout=5)

    def wait_forever(self) -> None:
        while self.thread.is_alive():
            time.sleep(0.2)



def make_handler(
    project_root: str | Path,
    config: DashboardConfig | None = None,
    template_html: str | None = None,
    runtime: DashboardRuntime | None = None,
) -> type[BaseHTTPRequestHandler]:
    root = Path(project_root).resolve()
    app_runtime = runtime or DashboardRuntime(
        workspace_root=root,
        launch_root=root.parent,
        active_root=root,
    )
    cfg = config or DashboardConfig()
    html = template_html if template_html is not None else _load_template()
    refresh_seconds = min(60.0, max(1.0, float(cfg.refresh_seconds or 3.0)))
    html = html.replace("__FASTMDX_REFRESH_SECONDS__", f"{refresh_seconds:g}")

    class LiveDashboardHandler(BaseHTTPRequestHandler):
        server_version = "FastMDXLive/1.0"

        def do_GET(self) -> None:  # noqa: N802 - stdlib API
            try:
                self._dispatch()
            except Exception as exc:  # noqa: BLE001 — dashboard must never crash the sim
                logger.warning("dashboard route failed: %s", exc)
                self.send_error(500, "Dashboard internal error")

        def do_POST(self) -> None:  # noqa: N802 - stdlib API
            try:
                self._dispatch_post()
            except Exception as exc:  # noqa: BLE001 — launcher errors stay local
                logger.warning("dashboard POST route failed: %s", exc)
                self._send_json(
                    {"ok": False, "error": "Dashboard request failed."},
                    status=500,
                )

        def _dispatch(self) -> None:
            parsed = urlparse(self.path)
            path = parsed.path
            root = app_runtime.data_root()
            if path in {"/", "/index", "/results", "/live"}:
                self._send_html(html)
                return
            if path == "/api/app-state" or path == "/api/launcher/state":
                self._send_json(app_runtime.snapshot())
                return
            if path == "/api/launcher/defaults":
                self._send_json(launcher_defaults())
                return
            if path == "/api/status":
                status = read_status(root)
                metrics = read_metrics(root)
                payload = {
                    "status": status,
                    "health": analyze_health(status, metrics),
                }
                self._send_json(payload)
                return
            if path == "/api/metrics":
                self._send_json({"metrics": read_metrics(root)})
                return
            if path == "/api/events":
                self._send_json({"events": read_events(root)})
                return
            if path == "/api/artifacts" or path == "/api/files":
                self._send_json({"artifacts": _artifact_records(root)})
                return
            if path == "/api/results" or path == "/api/analyses":
                self._send_json(_results_payload(root))
                return
            if path == "/api/protein-preview":
                regenerate = parsed.query == "regenerate=1"
                self._send_json(protein_preview_payload(root, regenerate=regenerate))
                return
            if path == "/api/structure-info":
                self._send_json(_structure_info_payload(root, cfg))
                return
            if path == "/api/ligands":
                self._send_json(_ligands_payload(root, cfg))
                return
            if path == "/api/live-frame-index":
                sim_dir = root / "simulation"
                self._send_json(read_live_frame_index(sim_dir))
                return
            if path == "/api/live-coordinates":
                self._send_json(_live_coordinates_payload(root))
                return
            if path == "/api/playback-info":
                max_frames = cfg.max_browser_frames
                if "max=" in parsed.query:
                    try:
                        max_frames = int(parsed.query.split("max=", 1)[1].split("&")[0])
                    except ValueError:
                        pass
                force = "force=1" in parsed.query
                sim_manifest = _load_json(root / "simulation" / "simulation_parameters.json")
                duration_ns = sim_manifest.get("duration_ns_actual")
                self._send_json(playback_info(
                    root,
                    max_browser_frames=max_frames,
                    simulation_time_ns_total=duration_ns,
                    force=force,
                ))
                return
            if path == "/api/open-output":
                opened, detail = _open_local_path(root)
                self._send_json({
                    "opened": opened,
                    "path": str(root),
                    "detail": detail,
                })
                return
            if path == "/analysis-figures-svg.zip":
                self._send_svg_bundle(root)
                return
            if path == "/structure/topology.pdb":
                self._send_structure(root)
                return
            if path == "/structure/live-frame.pdb":
                self._send_live_frame(root)
                return
            if path == "/structure/playback.pdb":
                self._send_playback(root)
                return
            if path.startswith("/static/"):
                self._send_static_asset(path.removeprefix("/static/"))
                return
            if path.startswith("/artifacts/"):
                self._send_artifact(
                    root,
                    path.removeprefix("/artifacts/"),
                    download="download=1" in parsed.query,
                )
                return
            self.send_error(404, "Not found")

        def _dispatch_post(self) -> None:
            parsed = urlparse(self.path)
            path = parsed.path
            payload = self._read_json_body()
            if path == "/api/launcher/validate":
                result = validate_launcher_payload(payload)
                if result.get("valid"):
                    config_payload = result["config"]
                    output_dir = app_runtime.launch_root / str(config_payload["run_name"])
                    result["output"] = str(output_dir)
                    result["command"] = build_launcher_command(config_payload, output_dir)
                self._send_json(result, status=200 if result.get("valid") else 422)
                return
            if path == "/api/launcher/launch":
                dashboard_url = self.headers.get("Origin")
                result = app_runtime.launch(payload, dashboard_url=dashboard_url)
                self._send_json(result, status=201 if result.get("launched") else 422)
                return
            if path == "/api/launcher/stop":
                self._send_json(app_runtime.stop())
                return
            self.send_error(404, "Not found")

        def log_message(self, format: str, *args: Any) -> None:  # noqa: A002
            return

        # ---- Generic response helpers ----
        def _read_json_body(self) -> dict[str, Any]:
            try:
                length = int(self.headers.get("Content-Length", "0") or 0)
            except ValueError:
                length = 0
            if length <= 0:
                return {}
            if length > 1_000_000:
                raise ValueError("Request body is too large")
            raw = self.rfile.read(length)
            data = json.loads(raw.decode("utf-8"))
            if not isinstance(data, dict):
                raise ValueError("JSON body must be an object")
            return data

        def _send_json(self, payload: dict[str, Any], *, status: int = 200) -> None:
            body = json.dumps(payload, default=str).encode("utf-8")
            self.send_response(status)
            self.send_header("Content-Type", "application/json; charset=utf-8")
            self.send_header("Cache-Control", "no-store")
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            self.wfile.write(body)

        def _send_html(self, html_text: str) -> None:
            body = html_text.encode("utf-8")
            self.send_response(200)
            self.send_header("Content-Type", "text/html; charset=utf-8")
            self.send_header("Cache-Control", "no-store")
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            self.wfile.write(body)

        def _send_artifact(
            self,
            root: Path,
            raw_rel: str,
            *,
            download: bool = False,
        ) -> None:
            try:
                target = (root / unquote(raw_rel)).resolve()
                target.relative_to(root)
            except ValueError:
                self.send_error(403, "Artifact path is outside the output directory")
                return
            if not target.is_file():
                self.send_error(404, "Artifact not found")
                return
            content_type = mimetypes.guess_type(target.name)[0] or "application/octet-stream"
            data = target.read_bytes()
            self.send_response(200)
            self.send_header("Content-Type", content_type)
            self.send_header("Cache-Control", "no-store")
            if download:
                safe_name = target.name.replace('"', "")
                self.send_header(
                    "Content-Disposition",
                    f'attachment; filename="{safe_name}"',
                )
            self.send_header("Content-Length", str(len(data)))
            self.end_headers()
            self.wfile.write(data)

        def _send_svg_bundle(self, root: Path) -> None:
            """Download every generated SVG analysis/report figure as one ZIP."""
            svg_paths = _svg_figure_paths(root)
            if not svg_paths:
                self.send_error(404, "No SVG figures are available yet")
                return
            buffer = io.BytesIO()
            with zipfile.ZipFile(buffer, mode="w", compression=zipfile.ZIP_DEFLATED) as archive:
                for source in svg_paths:
                    try:
                        relative = source.relative_to(root).as_posix()
                        archive.write(source, arcname=relative)
                    except (OSError, ValueError):
                        continue
            data = buffer.getvalue()
            self.send_response(200)
            self.send_header("Content-Type", "application/zip")
            self.send_header("Cache-Control", "no-store")
            self.send_header(
                "Content-Disposition",
                'attachment; filename="fastmdxplora_svg_figures.zip"',
            )
            self.send_header("Content-Length", str(len(data)))
            self.end_headers()
            self.wfile.write(data)

        def _send_structure(self, root: Path) -> None:
            target = find_structure(root)
            if target is None or not target.is_file():
                self.send_error(404, "Structure not found")
                return
            data = target.read_bytes()
            self.send_response(200)
            self.send_header("Content-Type", "chemical/x-pdb; charset=utf-8")
            self.send_header("Cache-Control", "no-store")
            self.send_header("Content-Length", str(len(data)))
            self.end_headers()
            self.wfile.write(data)

        def _send_live_frame(self, root: Path) -> None:
            sim_dir = root / "simulation"
            live_path = sim_dir / "live_frame.pdb"
            if not live_path.is_file():
                self.send_error(404, "Live frame not available")
                return
            try:
                data = live_path.read_bytes()
            except OSError:
                self.send_error(404, "Live frame not available")
                return
            self.send_response(200)
            self.send_header("Content-Type", "chemical/x-pdb; charset=utf-8")
            self.send_header("Cache-Control", "no-store")
            self.send_header("Content-Length", str(len(data)))
            self.end_headers()
            self.wfile.write(data)

        def _send_playback(self, root: Path) -> None:
            sim_dir = root / "simulation"
            playback_path = sim_dir / "playback.pdb"
            if not playback_path.is_file():
                # Auto-generate on demand so the dashboard never needs
                # to know whether the simulation phase has finished.
                playback_info(root, max_browser_frames=cfg.max_browser_frames)
            if not playback_path.is_file():
                self.send_error(404, "Playback not available")
                return
            try:
                data = playback_path.read_bytes()
            except OSError:
                self.send_error(404, "Playback not available")
                return
            self.send_response(200)
            self.send_header("Content-Type", "chemical/x-pdb; charset=utf-8")
            self.send_header("Cache-Control", "no-store")
            self.send_header("Content-Length", str(len(data)))
            self.end_headers()
            self.wfile.write(data)

        def _send_static_asset(self, raw_name: str) -> None:
            static_root = Path(__file__).with_name("static")
            target = (static_root / unquote(raw_name)).resolve()
            try:
                target.relative_to(static_root.resolve())
            except ValueError:
                self.send_error(404, "Static asset not found")
                return
            if not target.is_file():
                self.send_error(404, "Static asset not found")
                return
            content_type = mimetypes.guess_type(target.name)[0] or "application/octet-stream"
            data = target.read_bytes()
            self.send_response(200)
            self.send_header("Content-Type", content_type)
            self.send_header("Cache-Control", "no-store")
            self.send_header("Content-Length", str(len(data)))
            self.end_headers()
            self.wfile.write(data)

    return LiveDashboardHandler


def _load_template() -> str:
    try:
        return DASHBOARD_TEMPLATE_PATH.read_text(encoding="utf-8")
    except OSError:
        return (
            "<!doctype html><html><head><meta charset=\"utf-8\"><title>"
            "FastMDXplora Live Dashboard</title></head><body>"
            "<h1>Dashboard template unavailable</h1>"
            "<p>The dashboard HTML template could not be loaded. "
            "Reinstall the fastmdxplora package or run from the editable checkout.</p>"
            "</body></html>"
        )


def serve_dashboard(
    *,
    output: str | Path,
    host: str = "127.0.0.1",
    port: int = 8765,
    config: DashboardConfig | None = None,
    home_mode: bool = False,
    launch_root: str | Path | None = None,
) -> None:
    session = start_dashboard_session(
        output=output,
        host=host,
        port=port,
        config=config,
        home_mode=home_mode,
        launch_root=launch_root,
    )
    print(f"Live dashboard running at {session.url}")
    if session.port_was_changed:
        print(
            f"Requested port {session.requested_port} was busy, "
            f"so FastMDXplora used {session.port}."
        )
    print(f"Watching: {session.root}")
    print("Press Ctrl+C to stop.")
    try:
        session.wait_forever()
    except KeyboardInterrupt:
        pass
    finally:
        session.stop()


def start_dashboard_session(
    *,
    output: str | Path,
    host: str = "127.0.0.1",
    port: int = 8765,
    max_port_tries: int = 50,
    config: DashboardConfig | None = None,
    home_mode: bool = False,
    launch_root: str | Path | None = None,
) -> DashboardSession:
    """Start the local dashboard server in a background thread."""
    root = Path(output).resolve()
    root.mkdir(parents=True, exist_ok=True)
    runtime = DashboardRuntime(
        workspace_root=root,
        launch_root=Path(launch_root).resolve() if launch_root is not None else root.parent,
        active_root=None if home_mode else root,
    )
    requested_port = int(port)
    candidates = [0] if requested_port == 0 else range(requested_port, requested_port + max_port_tries)
    last_error: OSError | None = None
    for candidate in candidates:
        try:
            handler = make_handler(root, config=config, runtime=runtime)
            server = ThreadingHTTPServer((host, int(candidate)), handler)
        except OSError as exc:
            last_error = exc
            continue
        actual_port = int(server.server_address[1])
        thread = threading.Thread(
            target=server.serve_forever,
            name=f"FastMDXLiveDashboard:{actual_port}",
            daemon=True,
        )
        thread.start()
        return DashboardSession(
            server=server,
            thread=thread,
            root=root,
            host=host,
            port=actual_port,
            requested_port=requested_port,
            config=config,
            runtime=runtime,
        )
    if last_error is not None:
        raise last_error
    raise OSError("No dashboard ports were available")


def start_test_server(
    project_root: str | Path,
    *,
    config: DashboardConfig | None = None,
    home_mode: bool = False,
) -> tuple[ThreadingHTTPServer, str]:
    root = Path(project_root).resolve()
    runtime = DashboardRuntime(
        workspace_root=root,
        launch_root=root.parent,
        active_root=None if home_mode else root,
    )
    handler = make_handler(root, config=config, runtime=runtime)
    server = ThreadingHTTPServer(("127.0.0.1", 0), handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    return server, f"http://127.0.0.1:{server.server_address[1]}"


def _artifact_records(root: Path) -> list[dict[str, str]]:
    records: list[dict[str, str]] = []
    if not root.exists():
        return records
    for path in sorted(p for p in root.rglob("*") if p.is_file()):
        try:
            rel = path.relative_to(root).as_posix()
        except ValueError:
            continue
        if "__pycache__" in path.parts:
            continue
        records.append(
            {
                "path": rel,
                "name": path.name,
                "href": f"/artifacts/{rel}?v={int(path.stat().st_mtime)}",
                "download_href": (
                    f"/artifacts/{rel}?download=1&v={int(path.stat().st_mtime)}"
                ),
                "size": str(path.stat().st_size),
                "mtime": str(path.stat().st_mtime),
                "display_path": _compact_path(rel),
                "absolute_path": str(path.resolve()),
            }
        )
    return records


def _compact_path(path: str, *, keep: int = 2) -> str:
    parts = Path(path).parts
    if len(parts) <= keep:
        return path
    return ".../" + "/".join(parts[-keep:])


def _results_payload(root: Path) -> dict[str, Any]:
    artifacts = _artifact_records(root)
    by_path = {record["path"]: record for record in artifacts}
    manifest = _load_json(root / "manifest.json")
    analysis_manifest = _load_json(root / "analysis" / "analysis_manifest.json")
    sim_manifest = _load_json(root / "simulation" / "simulation_parameters.json")
    plots = _plot_records(artifacts)
    report_suffixes = {".md", ".html", ".htm", ".pdf", ".pptx", ".zip"}
    reports = [
        record
        for record in artifacts
        if record["path"].startswith("report/")
        and len(Path(record["path"]).parts) == 2
        and Path(record["path"]).suffix.lower() in report_suffixes
    ]
    reports.sort(key=lambda record: _report_order(record["path"]))
    dashboard = by_path.get("report/dashboard.html")
    summary = _summary_records(root, manifest, analysis_manifest, sim_manifest)
    setup_manifest = _load_json(root / "setup" / "setup_parameters.json")
    system = _system_info(root, manifest, analysis_manifest, sim_manifest)
    return {
        "refreshed_at": _iso_now(),
        "has_analysis": any(record["path"].startswith("analysis/") for record in artifacts),
        "has_report": any(record["path"].startswith("report/") for record in artifacts),
        "output_dir": str(root),
        "run_title": _run_title(root, manifest),
        "summary": summary,
        "system": system,
        "setup": _setup_details(setup_manifest),
        "simulation": _simulation_details(sim_manifest, read_status(root)),
        "phases": _phase_records(manifest),
        "analyses": _analysis_records(analysis_manifest, artifacts, plots),
        "dashboard": dashboard,
        "plots": plots,
        "key_plots": [plot for plot in plots if plot["title"] in KEY_PLOT_TITLES][:6],
        "svg_figure_count": len(_svg_figure_paths(root)),
        "svg_bundle_href": "/analysis-figures-svg.zip",
        "reports": reports,
        "artifacts": artifacts,
    }


def _structure_info_payload(
    root: Path, config: DashboardConfig
) -> dict[str, Any]:
    structure_path = find_structure(root)
    if structure_path is None:
        return {
            "valid": False,
            "reason": "missing",
            "ligand_resnames": [],
            "ligand_instances": [],
        }
    info = count_structure(structure_path)
    if not info.get("valid"):
        return dict(
            info,
            structure_available=True,
            structure_url="/structure/topology.pdb",
            ligand_resnames=[],
            ligand_instances=[],
        )
    info["atoms_by_resname"] = ligand_atom_counts(structure_path)
    info["explicit_ligand"] = config.ligand_resname
    try:
        info["structure_path"] = structure_path.relative_to(root).as_posix()
    except ValueError:
        info["structure_path"] = structure_path.as_posix()
    info["structure_url"] = (
        f"/structure/topology.pdb?v={int(structure_path.stat().st_mtime_ns)}"
    )
    info["structure_available"] = True
    return info


def _ligands_payload(
    root: Path, config: DashboardConfig
) -> dict[str, Any]:
    structure_path = find_structure(root)
    if structure_path is None or not structure_path.is_file():
        return {"ligands": [], "explicit": normalise_ligand_resname(config.ligand_resname), "valid": False}
    info = count_structure(structure_path)
    instances = []
    if info.get("valid"):
        explicit = [config.ligand_resname] if config.ligand_resname else None
        detected = detect_ligands(
            # Reconstruct (chain, resname, resi) keys from info — we
            # already collected them in count_structure. To avoid a
            # second PDB walk we accept that callers that want full
            # ligand IDs receive them via /api/structure-info.
            (
                (str(ins.get("chain", "A")), str(ins.get("resname", "")), str(ins.get("resi", "")))
                for ins in info.get("ligand_instances", [])
                if isinstance(ins, dict)
            ),
            explicit=explicit,
            include_cofactors=config.include_cofactors,
        )
        instances = [
            {
                "chain": inst["chain"],
                "resname": inst["resname"],
                "resi": inst["resi"],
                "explicit": inst["explicit"],
            }
            for inst in detected["instances"]
        ]
    return {
        "ligands": instances,
        "resnames": sorted({inst["resname"] for inst in instances}),
        "explicit": normalise_ligand_resname(config.ligand_resname),
        "include_cofactors": config.include_cofactors,
        "valid": True,
    }


def _live_coordinates_payload(root: Path) -> dict[str, Any]:
    sim_dir = root / "simulation"
    index = read_live_frame_index(sim_dir)
    return {
        "live_frame_exists": live_frame_exists(sim_dir),
        "index": index,
        "available": bool(index.get("live_frame_available")),
    }


def _load_json(path: Path) -> dict[str, Any]:
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {}
    return data if isinstance(data, dict) else {}


def _summary_records(
    root: Path,
    manifest: dict[str, Any],
    analysis_manifest: dict[str, Any],
    sim_manifest: dict[str, Any],
) -> list[dict[str, str]]:
    phase_statuses = [
        str(item.get("status") or "unknown").lower()
        for item in manifest.get("phases", [])
        if isinstance(item, dict)
    ]
    if phase_statuses and all(value in {"ok", "completed", "skipped"} for value in phase_statuses):
        status = "completed"
    elif any(value in {"error", "failed"} for value in phase_statuses):
        status = "failed"
    elif phase_statuses:
        status = "in progress"
    else:
        status = "not available"
    frames = _first_present(
        analysis_manifest.get("n_frames"),
        sim_manifest.get("n_production_frames"),
    )
    atoms = analysis_manifest.get("n_atoms")
    sim_time = sim_manifest.get("duration_ns_actual")
    live_status = read_status(root)
    temperature = _first_present(
        live_status.get("target_temperature_K"),
        _last_metric_value(root, "temperature"),
    )
    latest = _first_present(live_status.get("status"), status)
    return [
        {"label": "Project status", "value": status},
        {"label": "Latest status", "value": str(latest or "not available")},
        {"label": "Frames", "value": _display_value(frames)},
        {"label": "Atoms", "value": _display_value(atoms)},
        {
            "label": "Simulation time",
            "value": f"{sim_time} ns" if sim_time is not None else "not available",
        },
        {
            "label": "Temperature",
            "value": f"{temperature} K" if temperature not in (None, "") else "not available",
        },
    ]


def _last_metric_value(root: Path, field: str) -> str | None:
    metrics = read_metrics(root, limit=1)
    if not metrics:
        return None
    value = metrics[-1].get(field)
    return str(value) if value not in (None, "") else None


def _system_info(
    root: Path,
    manifest: dict[str, Any],
    analysis_manifest: dict[str, Any],
    sim_manifest: dict[str, Any],
) -> dict[str, str]:
    live_status = read_status(root)
    params = sim_manifest.get("parameters")
    if not isinstance(params, dict):
        params = sim_manifest
    return {
        "system": _display_value(manifest.get("system")),
        "output_folder": root.as_posix(),
        "atoms": _display_value(analysis_manifest.get("n_atoms")),
        "frames": _display_value(
            _first_present(
                analysis_manifest.get("n_frames"),
                sim_manifest.get("n_production_frames"),
            )
        ),
        "platform": _display_value(
            _first_present(
                live_status.get("platform"),
                sim_manifest.get("platform_used"),
                sim_manifest.get("platform"),
                params.get("platform"),
            )
        ),
        "timestep_fs": _display_value(
            _first_present(live_status.get("timestep_fs"), params.get("timestep_fs"))
        ),
        "checkpoint": _display_value(
            _first_present(
                live_status.get("current_checkpoint_path"),
                live_status.get("checkpoint_path"),
            )
        ),
    }


def _run_title(root: Path, manifest: dict[str, Any]) -> str:
    options = manifest.get("options") if isinstance(manifest.get("options"), dict) else {}
    report_options = options.get("report") if isinstance(options.get("report"), dict) else {}
    title = (
        report_options.get("title")
        or report_options.get("report_title")
        or manifest.get("title")
    )
    if title:
        return str(title)
    system = manifest.get("system")
    if system:
        return f"{Path(str(system)).stem} · FastMDXplora"
    return root.name or "FastMDXplora Live"


def _setup_details(setup_manifest: dict[str, Any]) -> dict[str, Any]:
    params = setup_manifest.get("parameters")
    if not isinstance(params, dict):
        params = setup_manifest
    forcefield = setup_manifest.get("resolved_forcefield")
    if not isinstance(forcefield, dict):
        forcefield = {}
    ligand = forcefield.get("ligand") or params.get("ligand")
    if isinstance(ligand, dict):
        ligand = ligand.get("name") or ligand.get("files")
    return {
        "ph": _first_present(params.get("ph"), params.get("pH")),
        "ion_concentration_M": _first_present(
            params.get("ion_concentration_M"),
            params.get("ion_concentration_m"),
        ),
        "force_field": _first_present(
            forcefield.get("name"),
            forcefield.get("protein_forcefield"),
            forcefield.get("forcefield"),
            params.get("forcefield"),
            params.get("force_field"),
        ),
        "water_model": _first_present(
            forcefield.get("water_model"),
            params.get("water_model"),
        ),
        "ligand": ligand,
    }


def _simulation_details(
    sim_manifest: dict[str, Any],
    live_status: dict[str, Any],
) -> dict[str, Any]:
    params = sim_manifest.get("parameters")
    if not isinstance(params, dict):
        params = sim_manifest
    return {
        "platform": _first_present(
            live_status.get("platform"),
            sim_manifest.get("platform_used"),
            sim_manifest.get("platform"),
            params.get("platform"),
        ),
        "precision": _first_present(
            live_status.get("precision"),
            params.get("precision"),
        ),
        "temperature_K": _first_present(
            live_status.get("target_temperature_K"),
            params.get("temperature_K"),
        ),
        "timestep_fs": _first_present(
            live_status.get("timestep_fs"),
            params.get("timestep_fs"),
        ),
        "friction_per_ps": params.get("friction_per_ps"),
        "pressure_bar": _first_present(
            params.get("pressure_bar"),
            params.get("pressure_atm"),
        ),
        "duration_ns_actual": _first_present(
            sim_manifest.get("duration_ns_actual"),
            params.get("duration_ns"),
        ),
        "n_production_frames": _first_present(
            sim_manifest.get("n_production_frames"),
            live_status.get("current_frame_count"),
        ),
    }


def _analysis_records(
    analysis_manifest: dict[str, Any],
    artifacts: list[dict[str, str]],
    plots: list[dict[str, str]],
) -> list[dict[str, Any]]:
    results = analysis_manifest.get("results")
    if not isinstance(results, dict):
        if analysis_manifest.get("status"):
            return [
                {
                    "name": "analysis",
                    "title": "Analysis",
                    "status": str(analysis_manifest.get("status")),
                    "message": str(analysis_manifest.get("note") or ""),
                    "artifacts": [],
                    "plot": None,
                }
            ]
        return []

    records: list[dict[str, Any]] = []
    for name, raw in results.items():
        data = raw if isinstance(raw, dict) else {}
        prefix = f"analysis/{name}/"
        related = [record for record in artifacts if record["path"].startswith(prefix)]
        plot = next(
            (
                item
                for item in plots
                if item["path"].startswith(prefix)
                or _normalise_analysis_name(item["title"]) == _normalise_analysis_name(str(name))
            ),
            None,
        )
        records.append(
            {
                "name": str(name),
                "title": PLOT_TITLE_ALIASES.get(str(name).lower(), _humanize_stem(str(name))),
                "status": str(data.get("status") or "unknown"),
                "message": str(data.get("message") or ""),
                "started_at": data.get("started_at"),
                "finished_at": data.get("finished_at"),
                "artifacts": related,
                "plot": plot,
            }
        )
    return records


def _normalise_analysis_name(value: str) -> str:
    return "".join(ch for ch in value.lower() if ch.isalnum())


def _report_order(path: str) -> tuple[int, str]:
    preferred = {
        "report/dashboard.html": 0,
        "report/report.md": 1,
        "report/slides.pptx": 2,
        "report/project_bundle.zip": 3,
    }
    return preferred.get(path, 99), path


def _first_present(*values: Any) -> Any:
    return next((value for value in values if value not in (None, "")), None)


def _display_value(value: Any) -> str:
    return "not available" if value in (None, "") else str(value)


def _phase_records(manifest: dict[str, Any]) -> list[dict[str, str]]:
    seen: dict[str, str] = {}
    for phase in manifest.get("phases", []):
        if isinstance(phase, dict):
            name = str(phase.get("name") or "").lower()
            if name:
                seen[name] = str(phase.get("status") or "unknown")
    records = []
    for name in ("setup", "simulation", "analysis", "report"):
        status = seen.get(name, "not run")
        records.append({"name": name.title() if name != "simulation" else "Simulation", "status": status})
    return records


def _plot_records(artifacts: list[dict[str, str]]) -> list[dict[str, str]]:
    """Return one display record per scientific figure, paired with SVG.

    PNG is preferred for fast browser previews while a matching true-vector
    SVG is exposed for publication download.  Pairing by normalized title also
    handles report/dashboard_assets copies and per-analysis artifact names.
    """
    image_suffixes = {".png", ".jpg", ".jpeg", ".gif", ".svg"}
    candidates = [
        record for record in artifacts
        if Path(record["path"]).suffix.lower() in image_suffixes
        and (
            record["path"].startswith("report/dashboard_assets/")
            or record["path"].startswith("analysis/")
        )
    ]
    grouped: dict[str, list[dict[str, str]]] = {}
    for record in candidates:
        title = _plot_title(record["path"])
        if title:
            grouped.setdefault(title, []).append(record)

    results: list[dict[str, str]] = []
    for title, records in grouped.items():
        display = min(records, key=lambda item: _plot_record_priority(item["path"]))
        svg_candidates = [item for item in records if Path(item["path"]).suffix.lower() == ".svg"]
        svg = min(svg_candidates, key=lambda item: _plot_priority(item["path"])) if svg_candidates else None
        mode = "dashboard view" if display["path"].startswith("report/dashboard_assets/") else "artifact fallback"
        enriched = dict(display)
        enriched.update(
            {
                "title": title,
                "category": PLOT_CATEGORY_BY_TITLE.get(title, _plot_category(display["path"])),
                "mode": mode,
            }
        )
        if svg is not None:
            enriched.update(
                {
                    "svg_path": svg["path"],
                    "svg_href": svg["href"],
                    "svg_download_href": svg["download_href"],
                }
            )
        results.append(enriched)
    return sorted(results, key=lambda item: (_category_order(item["category"]), item["title"]))


def _plot_record_priority(path: str) -> tuple[int, int]:
    suffix = Path(path).suffix.lower()
    suffix_priority = {".png": 0, ".jpg": 1, ".jpeg": 1, ".gif": 2, ".svg": 3}
    return _plot_priority(path), suffix_priority.get(suffix, 9)


def _plot_priority(path: str) -> int:
    return 0 if path.startswith("report/dashboard_assets/") else 1


def _plot_title(path: str) -> str | None:
    stem = Path(path).stem
    if stem in DASHBOARD_ASSET_TITLES:
        return DASHBOARD_ASSET_TITLES[stem]
    if stem in PLOT_TITLE_ALIASES:
        return PLOT_TITLE_ALIASES[stem]
    lower = stem.lower()
    for key, title in PLOT_TITLE_ALIASES.items():
        if lower == key or lower.endswith(f"_{key}") or key in lower:
            return title
    return _humanize_stem(stem)


def _humanize_stem(stem: str) -> str:
    return stem.replace("_", " ").replace("-", " ").strip().title()


def _plot_category(path: str) -> str:
    parts = Path(path).parts
    if "cluster" in parts or "dashboard_assets" in parts and "cluster" in path:
        return "Clustering"
    if any(part in {"rmsd", "rmsf", "rg", "hbonds", "hbond"} for part in parts):
        return "Core Metrics"
    if any(part in {"sasa", "ss", "dimred", "dihedrals", "qvalue"} for part in parts):
        return "Additional Analysis"
    return "Other"


def _category_order(category: str) -> int:
    order = {"Core Metrics": 0, "Additional Analysis": 1, "Clustering": 2, "Other": 3}
    return order.get(category, 99)


def _svg_figure_paths(root: Path) -> list[Path]:
    """Return generated SVG figures, excluding branding/static assets."""
    paths: list[Path] = []
    for base in (root / "analysis", root / "report" / "dashboard_assets", root / "report"):
        if not base.is_dir():
            continue
        for path in base.rglob("*.svg"):
            if path.is_file() and path not in paths:
                paths.append(path)
    return sorted(paths)


def _iso_now() -> str:
    from datetime import datetime, timezone

    return datetime.now(timezone.utc).isoformat()


def _open_local_path(path: Path) -> tuple[bool, str]:
    """Open ``path`` in the host file manager on a best-effort basis."""

    try:
        if os.name == "nt":
            os.startfile(str(path))  # type: ignore[attr-defined]
        elif sys.platform == "darwin":
            subprocess.Popen(["open", str(path)])
        else:
            subprocess.Popen(["xdg-open", str(path)])
    except Exception as exc:  # noqa: BLE001 - dashboard helper only
        return False, str(exc)
    return True, "opened"
