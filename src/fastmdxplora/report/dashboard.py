"""Static HTML dashboard for FastMDXplora report outputs."""

from __future__ import annotations

import csv
import json
import os
from dataclasses import dataclass
from datetime import datetime, timezone
from html import escape
from pathlib import Path
from typing import TYPE_CHECKING, Any
from urllib.parse import quote

from fastmdxplora.report.context import load_phase_context
from fastmdxplora.utils.logging import get_logger

if TYPE_CHECKING:
    from fastmdxplora.orchestrator import FastMDXplora

logger = get_logger("report.dashboard")

SECTION_ORDER: tuple[str, ...] = (
    "Core Metrics",
    "Additional Analysis",
    "SASA",
    "Secondary Structure",
    "Dimensionality Reduction",
    "Clustering",
    "Region Highlights",
    "Apo/Holo Comparison",
    "Other",
)

SECTION_ANCHORS: dict[str, str] = {
    "Core Metrics": "core-metrics",
    "Additional Analysis": "additional-analysis",
    "SASA": "sasa-section",
    "Secondary Structure": "secondary-structure-section",
    "Dimensionality Reduction": "dimensionality-reduction",
    "Clustering": "clustering-section",
    "Region Highlights": "region-highlights",
    "Apo/Holo Comparison": "apo-holo-comparison",
    "Other": "other-analysis",
}

ANALYSIS_SECTION_BY_FOLDER: dict[str, str] = {
    "rmsd": "Core Metrics",
    "rmsf": "Core Metrics",
    "rg": "Core Metrics",
    "hbonds": "Core Metrics",
    "sasa": "SASA",
    "ss": "Secondary Structure",
    "dimred": "Dimensionality Reduction",
    "cluster": "Clustering",
    "apo_holo": "Apo/Holo Comparison",
    "dihedrals": "Other",
    "qvalue": "Other",
}

DASHBOARD_ASSET_TITLE_ALIASES: dict[str, tuple[str, ...]] = {
    "RMSD": ("RMSD",),
    "RMSF": ("RMSF",),
    "Radius of gyration": ("Radius of gyration", "Rg"),
    "Hydrogen bonds": ("Hydrogen bonds", "H-bonds"),
    "Total SASA": ("SASA", "Total SASA"),
    "PCA": ("PCA", "Dimensionality reduction PCA"),
    "MDS": ("MDS", "Dimensionality reduction MDS"),
    "t-SNE": ("t-SNE", "Dimensionality reduction t-SNE"),
    "KMeans trajectory scatter": ("KMeans trajectory scatter", "Cluster KMeans"),
    "KMeans population plot": ("KMeans population plot", "KMeans populations"),
    "Hierarchical trajectory scatter": ("Hierarchical trajectory scatter",),
    "Hierarchical population plot": ("Hierarchical population plot",),
    "Hierarchical dendrogram": ("Hierarchical dendrogram",),
    "DBSCAN trajectory scatter": ("DBSCAN trajectory scatter",),
    "DBSCAN population plot": ("DBSCAN population plot",),
    "Secondary structure": ("Secondary structure", "SS heatmap"),
    "Fraction of native contacts": ("Fraction of native contacts", "Q-value"),
    "Dihedrals": ("Dihedrals",),
}

DASHBOARD_ASSET_SPECS: tuple[tuple[str, str, str, str, str, str, str], ...] = (
    (
        "RMSD",
        "analysis/rmsd/rmsd.dat",
        "rmsd_dashboard.png",
        "#35a7ff",
        "Frame",
        "RMSD (nm)",
        "line",
    ),
    (
        "RMSF",
        "analysis/rmsf/rmsf.dat",
        "rmsf_dashboard.png",
        "#d957c8",
        "Residue",
        "RMSF (nm)",
        "line",
    ),
    (
        "Radius of gyration",
        "analysis/rg/rg.dat",
        "rg_dashboard.png",
        "#57c45d",
        "Frame",
        "Rg (nm)",
        "line",
    ),
    (
        "Hydrogen bonds",
        "analysis/hbonds/hbonds.dat",
        "hbonds_dashboard.png",
        "#2ed3e6",
        "Frame",
        "H-bonds",
        "line",
    ),
    (
        "Total SASA",
        "analysis/sasa/sasa.dat",
        "sasa_dashboard.png",
        "#43c7b7",
        "Frame",
        "SASA (nm^2)",
        "line",
    ),
    (
        "PCA",
        "analysis/dimred/dimred_pca.dat",
        "pca_dashboard.png",
        "#35a7ff",
        "PC1",
        "PC2",
        "scatter",
    ),
    (
        "MDS",
        "analysis/dimred/dimred_mds.dat",
        "mds_dashboard.png",
        "#7cc66a",
        "MDS 1",
        "MDS 2",
        "scatter",
    ),
    (
        "t-SNE",
        "analysis/dimred/dimred_tsne.dat",
        "tsne_dashboard.png",
        "#d957c8",
        "t-SNE 1",
        "t-SNE 2",
        "scatter",
    ),
    (
        "KMeans trajectory scatter",
        "analysis/cluster/cluster_kmeans.dat",
        "kmeans_trajectory_dashboard.png",
        "#7cc66a",
        "Frame",
        "Cluster",
        "cluster",
    ),
    (
        "KMeans population plot",
        "analysis/cluster/cluster_kmeans.dat",
        "kmeans_population_dashboard.png",
        "#7cc66a",
        "Cluster",
        "Frames",
        "cluster_counts",
    ),
    (
        "Hierarchical trajectory scatter",
        "analysis/cluster/cluster_hierarchical.dat",
        "hierarchical_trajectory_dashboard.png",
        "#efb35e",
        "Frame",
        "Cluster",
        "cluster",
    ),
    (
        "Hierarchical population plot",
        "analysis/cluster/cluster_hierarchical.dat",
        "hierarchical_population_dashboard.png",
        "#efb35e",
        "Cluster",
        "Frames",
        "cluster_counts",
    ),
    (
        "Hierarchical dendrogram",
        "analysis/cluster/hierarchical_linkage.npy",
        "hierarchical_dendrogram_dashboard.png",
        "#efb35e",
        "Frame",
        "Distance",
        "dendrogram",
    ),
    (
        "DBSCAN trajectory scatter",
        "analysis/cluster/cluster_dbscan.dat",
        "dbscan_trajectory_dashboard.png",
        "#39b7c9",
        "Frame",
        "Cluster",
        "cluster",
    ),
    (
        "DBSCAN population plot",
        "analysis/cluster/cluster_dbscan.dat",
        "dbscan_population_dashboard.png",
        "#39b7c9",
        "Cluster",
        "Frames",
        "cluster_counts",
    ),
    (
        "Secondary structure",
        "analysis/ss/ss.dat",
        "ss_dashboard.png",
        "#35a7ff",
        "Residue",
        "Frame",
        "ss",
    ),
    (
        "Fraction of native contacts",
        "analysis/qvalue/qvalue.dat",
        "qvalue_dashboard.png",
        "#57c45d",
        "Frame",
        "Q",
        "line",
    ),
    (
        "Dihedrals",
        "analysis/dihedrals/dihedrals.dat",
        "dihedrals_dashboard.png",
        "#d957c8",
        "Phi (degrees)",
        "Psi (degrees)",
        "dihedrals",
    ),
)


@dataclass(frozen=True)
class DashboardCard:
    label: str
    value: str
    detail: str = ""
    kind: str = "neutral"


@dataclass(frozen=True)
class DashboardPanel:
    title: str
    source: str
    href: str
    original_source: str
    original_href: str
    mode: str
    summary: str = ""
    category: str = ""


@dataclass(frozen=True)
class DashboardLink:
    label: str
    href: str
    detail: str = ""


@dataclass(frozen=True)
class DashboardAsset:
    rel_path: str
    summary: str


@dataclass(frozen=True)
class DashboardSection:
    title: str
    anchor: str
    panels: list[DashboardPanel]


@dataclass(frozen=True)
class PhaseRow:
    name: str
    status: str
    detail: str


@dataclass(frozen=True)
class MetricRow:
    metric: str
    average: str
    stddev: str
    unit: str


def build_dashboard(
    *,
    orchestrator: "FastMDXplora",
    output_dir: Path,
    title: str,
    include_bundle_link: bool = False,
) -> list[str]:
    """Write ``dashboard.html`` and return artifact paths relative to output_dir."""
    project_root = orchestrator.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    manifest = _load_json(project_root / "manifest.json")
    analysis_manifest = _load_json(project_root / "analysis" / "analysis_manifest.json")
    sim_manifest = _load_json(project_root / "simulation" / "simulation_parameters.json")
    phase_context = load_phase_context(project_root)
    generated_at = datetime.now(timezone.utc)

    cards = _summary_cards(
        project_root=project_root,
        manifest=manifest,
        analysis_manifest=analysis_manifest,
        sim_manifest=sim_manifest,
    )
    dashboard_assets = _build_dashboard_assets(project_root, output_dir)
    sections = _analysis_sections(project_root, output_dir, dashboard_assets)
    links = _artifact_links(
        project_root,
        output_dir,
        sections=sections,
        include_bundle_link=include_bundle_link,
    )
    phase_rows = _phase_rows(manifest)
    metrics = _metric_rows(project_root, analysis_manifest)
    status = _project_status(manifest)
    live_html = _render_static_live_panel(project_root)
    phase_notice = ""
    if phase_context.is_analysis_from_existing_trajectory:
        phase_notice = (
            "Analysis/report workflow from existing trajectory. Setup and "
            "simulation were not run in this workflow."
        )

    html = _render_dashboard(
        title=title,
        system=str(manifest.get("system") or getattr(orchestrator, "system", "")),
        status=status,
        generated_at=generated_at,
        phase_notice=phase_notice,
        cards=cards,
        sections=sections,
        links=links,
        phase_rows=phase_rows,
        metrics=metrics,
        output_folder=project_root.as_posix(),
        live_html=live_html,
    )

    dashboard_path = output_dir / "dashboard.html"
    dashboard_path.write_text(html, encoding="utf-8")
    logger.debug("dashboard: wrote %s", dashboard_path)
    return ["dashboard.html"]


def _load_json(path: Path) -> dict[str, Any]:
    try:
        with path.open(encoding="utf-8") as fh:
            data = json.load(fh)
        return data if isinstance(data, dict) else {}
    except (FileNotFoundError, OSError, json.JSONDecodeError):
        return {}


def _summary_cards(
    *,
    project_root: Path,
    manifest: dict[str, Any],
    analysis_manifest: dict[str, Any],
    sim_manifest: dict[str, Any],
) -> list[DashboardCard]:
    phases = [p for p in manifest.get("phases", []) if isinstance(p, dict)]
    completed = [
        str(p.get("name"))
        for p in phases
        if p.get("status") == "ok" and p.get("name")
    ]
    status = _project_status(manifest)
    cards = [
        DashboardCard(
            "Project status",
            status.title(),
            "Recorded from manifest" if manifest else "No run manifest found",
            "good" if status == "ok" else "warn",
        ),
        DashboardCard(
            "Phases completed",
            str(len(completed)),
            ", ".join(completed) if completed else "not available",
        ),
    ]

    n_frames = analysis_manifest.get("n_frames")
    if n_frames is not None:
        cards.append(DashboardCard("Frames", _format_number(n_frames), "analysis metadata"))

    n_atoms = analysis_manifest.get("n_atoms")
    if n_atoms is not None:
        cards.append(DashboardCard("Atom count", _format_number(n_atoms), "topology metadata"))

    phase_context = load_phase_context(project_root)
    if phase_context.simulation_present:
        sim_params = sim_manifest.get("parameters", {})
        if isinstance(sim_params, dict):
            duration = sim_params.get("duration_ns")
            if duration is not None:
                cards.append(
                    DashboardCard(
                        "Simulation time",
                        f"{_format_number(duration)} ns",
                        "simulation parameters",
                    )
                )
        temperature = _average_temperature(project_root / "simulation" / "energy.csv")
        if temperature is not None:
            cards.append(
                DashboardCard(
                    "Temperature",
                    f"{temperature:.1f} K",
                    "energy log average",
                )
            )

    wall_time = _wall_time(phases)
    if wall_time:
        cards.append(DashboardCard("Wall time", wall_time, "recorded phase timestamps"))

    cards.append(
        DashboardCard("Output folder", project_root.as_posix(), "project root")
    )
    return cards


def _project_status(manifest: dict[str, Any]) -> str:
    phases = [p for p in manifest.get("phases", []) if isinstance(p, dict)]
    if not phases:
        return "unknown"
    if any(p.get("status") == "error" for p in phases):
        return "error"
    if all(p.get("status") in {"ok", "skipped"} for p in phases):
        return "ok"
    return "unknown"


def _phase_rows(manifest: dict[str, Any]) -> list[PhaseRow]:
    recorded = {
        str(p.get("name")): p
        for p in manifest.get("phases", [])
        if isinstance(p, dict) and p.get("name")
    }
    rows: list[PhaseRow] = []
    for name in ("setup", "simulation", "analysis", "report"):
        phase = recorded.get(name)
        if phase is None:
            rows.append(PhaseRow(name.title(), "not-run", "Not run"))
            continue
        raw_status = str(phase.get("status") or "unknown")
        if raw_status == "ok":
            detail = "Completed"
        elif raw_status == "error":
            detail = "Failed"
        elif raw_status == "skipped":
            detail = "Skipped"
        else:
            detail = raw_status.title()
        rows.append(PhaseRow(name.title(), raw_status, detail))
    return rows


def _metric_rows(project_root: Path, analysis_manifest: dict[str, Any]) -> list[MetricRow]:
    specs: tuple[tuple[str, str, str], ...] = (
        ("RMSD", "analysis/rmsd/rmsd.dat", "nm"),
        ("RMSF", "analysis/rmsf/rmsf.dat", "nm"),
        ("Radius of gyration", "analysis/rg/rg.dat", "nm"),
        ("H-bonds", "analysis/hbonds/hbonds.dat", "count"),
        ("SASA", "analysis/sasa/sasa.dat", "nm^2"),
    )
    rows: list[MetricRow] = []
    for label, rel, unit in specs:
        values = _numeric_series(project_root / rel)
        if not values:
            continue
        rows.append(
            MetricRow(
                metric=label,
                average=_format_metric_value(_mean(values)),
                stddev=_format_metric_value(_stddev(values)),
                unit=unit,
            )
        )

    n_frames = analysis_manifest.get("n_frames")
    if n_frames is not None:
        rows.append(MetricRow("Frame count", _format_number(n_frames), "not available", "frames"))
    n_atoms = analysis_manifest.get("n_atoms")
    if n_atoms is not None:
        rows.append(MetricRow("Atom count", _format_number(n_atoms), "not available", "atoms"))
    return rows


def _numeric_series(path: Path) -> list[float]:
    if not path.is_file():
        return []
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except OSError:
        return []
    rows: list[list[float]] = []
    for line in lines:
        text = line.strip()
        if not text or text.startswith("#"):
            continue
        parts = [part for part in text.replace(",", " ").split() if part]
        values: list[float] = []
        for part in parts:
            try:
                values.append(float(part))
            except ValueError:
                values = []
                break
        if values:
            rows.append(values)
    if not rows:
        return []
    if all(len(row) == 1 for row in rows):
        return [row[0] for row in rows]
    return [row[-1] for row in rows if row]


def _mean(values: list[float]) -> float:
    return sum(values) / len(values)


def _stddev(values: list[float]) -> float:
    if len(values) < 2:
        return 0.0
    avg = _mean(values)
    return (sum((value - avg) ** 2 for value in values) / len(values)) ** 0.5


def _format_metric_value(value: float) -> str:
    if abs(value) >= 100:
        return f"{value:,.1f}"
    if abs(value) >= 10:
        return f"{value:,.2f}"
    if abs(value) >= 1:
        return f"{value:,.3f}"
    return f"{value:,.4f}"


def _average_temperature(path: Path) -> float | None:
    if not path.is_file():
        return None
    values: list[float] = []
    try:
        with path.open(newline="", encoding="utf-8") as fh:
            reader = csv.DictReader(_strip_comment_prefix(fh))
            for row in reader:
                raw = row.get("Temperature (K)")
                if raw in (None, "", "--"):
                    continue
                try:
                    values.append(float(raw))
                except ValueError:
                    continue
    except OSError:
        return None
    if not values:
        return None
    return sum(values) / len(values)


def _strip_comment_prefix(lines):
    for line in lines:
        yield line[1:] if line.startswith("#") else line


def _wall_time(phases: list[dict[str, Any]]) -> str:
    starts: list[datetime] = []
    finishes: list[datetime] = []
    for phase in phases:
        start = _parse_datetime(phase.get("started_at"))
        finish = _parse_datetime(phase.get("finished_at"))
        if start:
            starts.append(start)
        if finish:
            finishes.append(finish)
    if not starts or not finishes:
        return ""
    seconds = max(0, int((max(finishes) - min(starts)).total_seconds()))
    if seconds < 60:
        return f"{seconds}s"
    minutes, rem = divmod(seconds, 60)
    if minutes < 60:
        return f"{minutes}m {rem}s"
    hours, minutes = divmod(minutes, 60)
    return f"{hours}h {minutes}m"


def _parse_datetime(value: object) -> datetime | None:
    if not value:
        return None
    try:
        return datetime.fromisoformat(str(value).replace("Z", "+00:00"))
    except ValueError:
        return None


def _build_dashboard_assets(
    project_root: Path,
    output_dir: Path,
) -> dict[str, DashboardAsset]:
    assets: dict[str, DashboardAsset] = {}
    asset_dir = output_dir / "dashboard_assets"
    if asset_dir.is_dir():
        for stale in asset_dir.glob("*_dashboard.png"):
            try:
                stale.unlink()
            except OSError:
                logger.warning("dashboard: could not remove stale asset %s", stale)
    for title, rel_data, filename, color, xlabel, ylabel, kind in DASHBOARD_ASSET_SPECS:
        data_path = project_root / rel_data
        if not data_path.is_file():
            continue
        asset_path = asset_dir / filename
        try:
            asset_dir.mkdir(parents=True, exist_ok=True)
            summary = _write_dashboard_chart(
                data_path=data_path,
                output_path=asset_path,
                kind=kind,
                color=color,
                xlabel=xlabel,
                ylabel=ylabel,
            )
        except Exception as exc:  # noqa: BLE001 -- fallback to original figure
            logger.warning("dashboard chart skipped for %s: %s", title, exc)
            continue
        assets[title] = DashboardAsset(
            rel_path=f"report/dashboard_assets/{filename}",
            summary=summary,
        )
    return assets


def _write_dashboard_chart(
    *,
    data_path: Path,
    output_path: Path,
    kind: str,
    color: str,
    xlabel: str,
    ylabel: str,
) -> str:
    import matplotlib

    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(5.6, 3.15), dpi=180)
    _style_dashboard_axes(fig, ax)

    if kind == "scatter":
        rows = _numeric_rows(data_path)
        if not rows or any(len(row) < 3 for row in rows):
            raise ValueError("PCA data must include frame and two components")
        frame = [row[0] for row in rows]
        x = [row[1] for row in rows]
        y = [row[2] for row in rows]
        points = ax.scatter(
            x,
            y,
            c=frame,
            cmap="viridis",
            s=28,
            edgecolors="none",
            alpha=0.95,
        )
        cbar = fig.colorbar(points, ax=ax, pad=0.02, fraction=0.05)
        cbar.ax.tick_params(colors="black", labelsize=7)
        cbar.outline.set_edgecolor("#666666")
        cbar.set_label("Frame", color="black", fontsize=8)
        summary = f"{len(rows)} frames"
    elif kind == "cluster":
        rows = _numeric_rows(data_path)
        if not rows or any(len(row) < 2 for row in rows):
            raise ValueError("cluster data must include frame and cluster columns")
        x = [row[0] for row in rows]
        y = [row[1] for row in rows]
        ax.step(x, y, where="mid", color=color, linewidth=1.8)
        ax.scatter(x, y, color="#39b7c9", s=18, edgecolors="none", alpha=0.9)
        clusters = sorted({int(value) for value in y})
        summary = f"{len(clusters)} clusters"
    elif kind == "cluster_counts":
        rows = _numeric_rows(data_path)
        if not rows or any(len(row) < 2 for row in rows):
            raise ValueError("cluster data must include frame and cluster columns")
        counts: dict[int, int] = {}
        for row in rows:
            cluster = int(row[1])
            counts[cluster] = counts.get(cluster, 0) + 1
        clusters = sorted(counts)
        values = [counts[cluster] for cluster in clusters]
        ax.bar(
            [str(cluster) for cluster in clusters],
            values,
            color=color,
            edgecolor="#b8f1ff",
            linewidth=0.4,
            alpha=0.88,
        )
        summary = f"{len(clusters)} clusters"
    elif kind == "dihedrals":
        rows = _numeric_rows(data_path)
        if not rows or any(len(row) < 4 for row in rows):
            raise ValueError("dihedral data must include frame, residue, phi, and psi")
        frames = [row[0] for row in rows]
        phi = [row[2] for row in rows]
        psi = [row[3] for row in rows]
        points = ax.scatter(
            phi,
            psi,
            c=frames,
            cmap="plasma",
            s=12,
            edgecolors="none",
            alpha=0.82,
        )
        ax.set_xlim(-180, 180)
        ax.set_ylim(-180, 180)
        cbar = fig.colorbar(points, ax=ax, pad=0.02, fraction=0.05)
        cbar.ax.tick_params(colors="black", labelsize=7)
        cbar.outline.set_edgecolor("#666666")
        cbar.set_label("Frame", color="black", fontsize=8)
        summary = f"{len(rows)} angles"
    elif kind == "ss":
        matrix, residues, frames = _secondary_structure_matrix(data_path)
        if not matrix:
            raise ValueError("secondary structure data is empty")
        from matplotlib.colors import ListedColormap

        cmap = ListedColormap(["#1b3552", "#35a7ff", "#57c45d", "#d957c8"])
        ax.imshow(matrix, aspect="auto", interpolation="nearest", cmap=cmap)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_xticks([0, max(0, len(residues) - 1)])
        ax.set_xticklabels([residues[0], residues[-1]])
        ax.set_yticks([0, max(0, len(frames) - 1)])
        ax.set_yticklabels([frames[0], frames[-1]])
        summary = f"{len(frames)} frames, {len(residues)} residues"
        _finish_dashboard_chart(fig, ax, output_path)
        return summary
    elif kind == "dendrogram":
        summary = _plot_dashboard_dendrogram(
            ax=ax,
            linkage_path=data_path,
            color=color,
            xlabel=xlabel,
            ylabel=ylabel,
        )
        _finish_dashboard_chart(fig, ax, output_path)
        return summary
    else:
        rows = _numeric_rows(data_path)
        if not rows:
            raise ValueError("numeric data is empty")
        if all(len(row) == 1 for row in rows):
            x = list(range(len(rows)))
            y = [row[0] for row in rows]
        else:
            x = [row[0] for row in rows]
            y = [row[-1] for row in rows]
        ax.plot(x, y, color=color, linewidth=1.8)
        ax.fill_between(x, y, min(y), color=color, alpha=0.10)
        summary = f"avg {_format_metric_value(_mean(y))}"

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    _finish_dashboard_chart(fig, ax, output_path)
    return summary


def _plot_dashboard_dendrogram(
    *,
    ax,
    linkage_path: Path,
    color: str,
    xlabel: str,
    ylabel: str,
) -> str:
    import numpy as np

    try:
        from scipy.cluster.hierarchy import dendrogram
    except ImportError as exc:  # pragma: no cover - optional dependency
        raise RuntimeError("SciPy is required for dashboard dendrogram rendering") from exc

    linkage_matrix = np.load(linkage_path)
    if linkage_matrix.ndim != 2 or linkage_matrix.shape[1] != 4:
        raise ValueError("hierarchical linkage data must be an n x 4 matrix")
    dendrogram(
        linkage_matrix,
        ax=ax,
        no_labels=True,
        color_threshold=None,
        link_color_func=lambda _node_id: color,
        above_threshold_color=color,
    )
    for collection in ax.collections:
        collection.set_color(color)
        collection.set_linewidth(1.2)
        collection.set_alpha(0.95)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_ylim(bottom=0)
    return f"{linkage_matrix.shape[0] + 1} frames"


def _style_dashboard_axes(fig, ax) -> None:
    """Apply publication-safe styling to report chart assets."""
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")
    for spine in ax.spines.values():
        spine.set_color("#333333")
    ax.tick_params(colors="black", labelsize=8)
    ax.xaxis.label.set_color("black")
    ax.yaxis.label.set_color("black")
    ax.xaxis.label.set_size(8)
    ax.yaxis.label.set_size(8)
    ax.grid(True, color="#d9d9d9", alpha=0.8, linewidth=0.6, linestyle="--")


def _finish_dashboard_chart(fig, ax, output_path: Path) -> None:
    fig.tight_layout(pad=0.75)
    fig.savefig(
        output_path,
        dpi=300,
        facecolor="white",
        edgecolor="white",
        transparent=False,
        bbox_inches="tight",
    )
    if output_path.suffix.lower() != ".svg":
        fig.savefig(
            output_path.with_suffix(".svg"),
            format="svg",
            facecolor="white",
            edgecolor="white",
            transparent=False,
            bbox_inches="tight",
        )
    import matplotlib.pyplot as plt

    plt.close(fig)


def _numeric_rows(path: Path) -> list[list[float]]:
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except OSError:
        return []
    rows: list[list[float]] = []
    for line in lines:
        text = line.strip()
        if not text or text.startswith("#"):
            continue
        values: list[float] = []
        for part in [part for part in text.replace(",", " ").split() if part]:
            try:
                values.append(float(part))
            except ValueError:
                values = []
                break
        if values:
            rows.append(values)
    return rows


def _secondary_structure_matrix(path: Path) -> tuple[list[list[int]], list[str], list[str]]:
    code_map = {"C": 0, "H": 1, "E": 2, "B": 2, "G": 1, "I": 1, "T": 3, "S": 3}
    try:
        with path.open(newline="", encoding="utf-8") as fh:
            reader = csv.reader(fh)
            header = next(reader)
            residues = header[1:]
            matrix: list[list[int]] = []
            frames: list[str] = []
            for row in reader:
                if len(row) < 2:
                    continue
                frames.append(row[0])
                matrix.append([code_map.get(value, 0) for value in row[1:]])
    except (OSError, StopIteration):
        return [], [], []
    return matrix, residues, frames


def _analysis_sections(
    project_root: Path,
    output_dir: Path,
    dashboard_assets: dict[str, DashboardAsset],
) -> list[DashboardSection]:
    grouped: dict[str, list[DashboardPanel]] = {section: [] for section in SECTION_ORDER}
    seen_sources: set[str] = set()

    for source in _discover_analysis_images(project_root):
        rel = source.relative_to(project_root).as_posix()
        if rel in seen_sources:
            continue
        seen_sources.add(rel)
        folder = source.relative_to(project_root).parts[1]
        section = ANALYSIS_SECTION_BY_FOLDER.get(folder, "Other")
        title = _figure_title_from_path(source)
        asset = _asset_for_title(title, dashboard_assets)
        if asset is not None:
            display_path = project_root / asset.rel_path
            display_rel = asset.rel_path
            href = _href(display_path, output_dir)
            mode = "dashboard view"
            summary = asset.summary
        else:
            display_rel = rel
            href = _href(source, output_dir)
            mode = "artifact fallback"
            summary = ""
        grouped[section].append(
            DashboardPanel(
                title=title,
                source=display_rel,
                href=href,
                original_source=rel,
                original_href=_href(source, output_dir),
                mode=mode,
                summary=summary,
                category=section,
            )
        )

    for source in _discover_report_images(project_root):
        rel = source.relative_to(project_root).as_posix()
        if rel in seen_sources:
            continue
        seen_sources.add(rel)
        title = _figure_title_from_path(source)
        section = (
            "Region Highlights"
            if "region" in source.stem.lower()
            else "Apo/Holo Comparison"
            if "apo" in source.stem.lower() or "holo" in source.stem.lower()
            else "Other"
        )
        grouped[section].append(
            DashboardPanel(
                title=title,
                source=rel,
                href=_href(source, output_dir),
                original_source=rel,
                original_href=_href(source, output_dir),
                mode="artifact fallback",
                summary="",
                category=section,
            )
        )

    grouped = _group_sparse_sections(grouped)
    sections: list[DashboardSection] = []
    for title in SECTION_ORDER:
        panels = sorted(grouped[title], key=lambda panel: _panel_sort_key(panel))
        if panels:
            sections.append(
                DashboardSection(
                    title=title,
                    anchor=SECTION_ANCHORS[title],
                    panels=panels,
                )
            )
    return sections


def _group_sparse_sections(
    grouped: dict[str, list[DashboardPanel]],
) -> dict[str, list[DashboardPanel]]:
    """Move small non-core/non-clustering sections into Additional Analysis."""
    merged = {section: list(panels) for section, panels in grouped.items()}
    additional: list[DashboardPanel] = list(merged.get("Additional Analysis", []))
    for section, panels in list(merged.items()):
        if section in {"Core Metrics", "Additional Analysis", "Clustering"}:
            continue
        if panels and len(panels) <= 3:
            additional.extend(panels)
            merged[section] = []
    merged["Additional Analysis"] = additional
    return merged


def _discover_analysis_images(project_root: Path) -> list[Path]:
    analysis_dir = project_root / "analysis"
    if not analysis_dir.is_dir():
        return []
    suffixes = {".png", ".jpg", ".jpeg", ".svg"}
    return sorted(
        path
        for path in analysis_dir.rglob("*")
        if path.is_file() and path.suffix.lower() in suffixes
    )


def _discover_report_images(project_root: Path) -> list[Path]:
    report_dir = project_root / "report"
    if not report_dir.is_dir():
        return []
    names = {
        "region_highlight_summary.png",
        "structure_region_highlights.png",
        "apo_holo_comparison.png",
    }
    return sorted(path for path in report_dir.iterdir() if path.name in names)


def _asset_for_title(
    title: str,
    dashboard_assets: dict[str, DashboardAsset],
) -> DashboardAsset | None:
    normalized = _normalize_label(title)
    for asset_title, aliases in DASHBOARD_ASSET_TITLE_ALIASES.items():
        if _normalize_label(asset_title) == normalized:
            return dashboard_assets.get(asset_title)
        if any(_normalize_label(alias) == normalized for alias in aliases):
            return dashboard_assets.get(asset_title)
    return None


def _figure_title_from_path(path: Path) -> str:
    stem = path.stem
    folder = path.parent.name
    special = {
        "rmsd": "RMSD",
        "rmsf": "RMSF",
        "rg": "Radius of gyration",
        "hbonds": "Hydrogen bonds",
        "ss": "Secondary structure",
        "sasa": "Total SASA",
        "sasa_heatmap": "Per-residue SASA heatmap",
        "sasa_by_residue": "Average per-residue SASA",
        "total_sasa": "Total SASA",
        "residue_sasa": "Per-residue SASA heatmap",
        "average_residue_sasa": "Average per-residue SASA",
        "dimred_pca": "PCA",
        "dimred_mds": "MDS",
        "dimred_tsne": "t-SNE",
        "cluster_kmeans": "KMeans trajectory scatter",
        "cluster_kmeans_counts": "KMeans population plot",
        "cluster_dbscan": "DBSCAN trajectory scatter",
        "cluster_dbscan_counts": "DBSCAN population plot",
        "cluster_hierarchical": "Hierarchical trajectory scatter",
        "cluster_hierarchical_counts": "Hierarchical population plot",
        "cluster_hierarchical_dendrogram": "Hierarchical dendrogram",
        "dbscan_pop": "DBSCAN population plot",
        "dbscan_traj_hist": "DBSCAN trajectory histogram",
        "dbscan_traj_scatter": "DBSCAN trajectory scatter",
        "dbscan_distance_matrix": "DBSCAN distance matrix",
        "kmeans_pop": "KMeans population plot",
        "kmeans_traj_hist": "KMeans trajectory histogram",
        "kmeans_traj_scatter": "KMeans trajectory scatter",
        "hierarchical_pop": "Hierarchical population plot",
        "hierarchical_traj_hist": "Hierarchical trajectory histogram",
        "hierarchical_traj_scatter": "Hierarchical trajectory scatter",
        "hierarchical_dendrogram": "Hierarchical dendrogram",
        "region_highlight_summary": "Region highlights",
        "structure_region_highlights": "Structure region highlights",
        "apo_holo_comparison": "Apo/Holo comparison",
        "qvalue": "Fraction of native contacts",
        "dihedrals": "Dihedrals",
    }
    if stem in special:
        return special[stem]
    if folder == "dimred" and stem.startswith("dimred_"):
        return _friendly_name(stem.replace("dimred_", ""))
    if folder == "cluster":
        return _friendly_name(stem)
    return _friendly_name(stem)


def _friendly_name(value: str) -> str:
    words = value.replace("-", "_").split("_")
    replacements = {
        "rmsd": "RMSD",
        "rmsf": "RMSF",
        "rg": "Rg",
        "sasa": "SASA",
        "pca": "PCA",
        "mds": "MDS",
        "tsne": "t-SNE",
        "dbscan": "DBSCAN",
        "kmeans": "KMeans",
        "ss": "SS",
    }
    return " ".join(replacements.get(word.lower(), word.title()) for word in words)


def _panel_sort_key(panel: DashboardPanel) -> tuple[int, str]:
    source = panel.original_source
    priorities = (
        "rmsd.png",
        "rmsf.png",
        "rg.png",
        "hbonds.png",
        "sasa.png",
        "total_sasa.png",
        "sasa_heatmap.png",
        "residue_sasa.png",
        "sasa_by_residue.png",
        "average_residue_sasa.png",
        "ss.png",
        "dimred_pca.png",
        "dimred_mds.png",
        "dimred_tsne.png",
        "dbscan",
        "kmeans",
        "hierarchical",
    )
    for index, pattern in enumerate(priorities):
        if pattern in source:
            return index, source
    return len(priorities), source


def _normalize_label(value: str) -> str:
    return "".join(char.lower() for char in value if char.isalnum())


def _artifact_links(
    project_root: Path,
    output_dir: Path,
    *,
    sections: list[DashboardSection],
    include_bundle_link: bool,
) -> list[DashboardLink]:
    candidates: list[tuple[str, str, str]] = [
        ("Dashboard HTML", "report/dashboard.html", "this file"),
        ("Markdown report", "report/report.md", "written report"),
        ("Slide deck", "report/slides.pptx", "presentation"),
        ("Project bundle", "report/project_bundle.zip", "shareable archive"),
        ("Run manifest", "manifest.json", "project provenance"),
        ("Analysis manifest", "analysis/analysis_manifest.json", "analysis provenance"),
        ("Analysis summary", "report/analysis_summary.png", "combined figure"),
        ("Analysis summary manifest", "report/analysis_summary_manifest.json", "figure provenance"),
        ("Region highlight manifest", "report/region_highlight_manifest.json", "region provenance"),
        ("Region highlight figure", "report/region_highlight_summary.png", "annotated RMSF"),
    ]
    for section in sections:
        for panel in section.panels:
            candidates.append((panel.title, panel.original_source, section.title))

    links: list[DashboardLink] = []
    seen: set[str] = set()
    for label, rel, detail in candidates:
        path = project_root / rel
        if rel in seen:
            continue
        future_current_run_artifact = rel == "report/dashboard.html" or (
            include_bundle_link and rel == "report/project_bundle.zip"
        )
        if not path.is_file() and not future_current_run_artifact:
            continue
        seen.add(rel)
        links.append(DashboardLink(label, _href(path, output_dir), detail))
    return links


def _href(path: Path, output_dir: Path) -> str:
    rel = os.path.relpath(path, output_dir).replace(os.sep, "/")
    return quote(rel, safe="/._-#")


def _format_number(value: object) -> str:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return str(value)
    if number.is_integer():
        return f"{int(number):,}"
    return f"{number:,.3g}"


def _render_dashboard(
    *,
    title: str,
    system: str,
    status: str,
    generated_at: datetime,
    phase_notice: str,
    cards: list[DashboardCard],
    sections: list[DashboardSection],
    links: list[DashboardLink],
    phase_rows: list[PhaseRow],
    metrics: list[MetricRow],
    output_folder: str,
    live_html: str,
) -> str:
    nav_html = _render_sidebar(sections, links)
    card_html = "\n".join(_render_card(card) for card in cards)
    sections_html = "\n".join(_render_section(section) for section in sections)
    if not sections_html:
        sections_html = (
            '<section class="empty-state panel-block">'
            "<h2>Analysis panels</h2>"
            "<p>No analysis figures are available in this run output.</p>"
            "</section>"
        )
    link_html = _render_output_list(links)
    quick_html = "\n".join(_render_quick_action(link) for link in _quick_action_links(links))
    phase_html = "\n".join(_render_phase_row(row) for row in phase_rows)
    metrics_html = "\n".join(_render_metric_row(row) for row in metrics)
    if not metrics_html:
        metrics_html = (
            '<tr><td colspan="4" class="table-empty">not available</td></tr>'
        )
    notice_html = (
        f'<div class="notice"><strong>Existing trajectory analysis</strong>'
        f"{escape(phase_notice)}</div>" if phase_notice else ""
    )
    generated = generated_at.strftime("%Y-%m-%d %H:%M UTC")

    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{escape(title)} - FastMDXplora Dashboard</title>
  <style>
    :root {{
      color-scheme: dark;
      --bg: #07101b;
      --sidebar: #071321;
      --panel: #101a28;
      --panel-2: #152437;
      --panel-3: #0c1724;
      --line: #22364b;
      --text: #edf4fb;
      --muted: #a7b5c6;
      --accent: #39b7c9;
      --accent-blue: #4d9df7;
      --accent-2: #7cc66a;
      --warn: #efb35e;
      --danger: #e35d6a;
      --shadow: rgba(0, 0, 0, 0.25);
    }}
    * {{ box-sizing: border-box; }}
    body {{
      margin: 0;
      background:
        radial-gradient(circle at top left, rgba(57, 183, 201, 0.12), transparent 34rem),
        linear-gradient(180deg, #07101b 0%, #0a1724 100%);
      color: var(--text);
      font-family: Inter, ui-sans-serif, system-ui, -apple-system, BlinkMacSystemFont,
        "Segoe UI", sans-serif;
      line-height: 1.5;
    }}
    a {{ color: inherit; }}
    .layout {{
      display: grid;
      grid-template-columns: 260px minmax(0, 1fr);
      min-height: 100vh;
    }}
    .sidebar {{
      position: sticky;
      top: 0;
      height: 100vh;
      overflow-y: auto;
      border-right: 1px solid var(--line);
      background: linear-gradient(180deg, rgba(7, 19, 33, 0.98), rgba(5, 13, 23, 0.98));
      padding: 22px 16px;
    }}
    .logo {{
      display: grid;
      grid-template-columns: 42px 1fr;
      gap: 12px;
      align-items: center;
      padding-bottom: 22px;
      border-bottom: 1px solid var(--line);
      margin-bottom: 20px;
    }}
    .mark {{
      display: grid;
      place-items: center;
      width: 42px;
      height: 42px;
      border-radius: 12px;
      border: 1px solid rgba(77, 157, 247, 0.5);
      background: rgba(77, 157, 247, 0.12);
      color: #7ed5ff;
      font-weight: 800;
    }}
    .logo-title {{
      font-size: 1.08rem;
      font-weight: 800;
      line-height: 1.1;
    }}
    .logo-subtitle {{
      color: var(--muted);
      font-size: 0.72rem;
      margin-top: 3px;
    }}
    .nav-section {{ margin: 18px 0; }}
    .nav-heading {{
      color: #7d8da2;
      font-size: 0.72rem;
      font-weight: 800;
      letter-spacing: 0.07em;
      text-transform: uppercase;
      margin: 0 0 8px;
    }}
    .nav-link {{
      display: flex;
      align-items: center;
      gap: 9px;
      min-height: 34px;
      padding: 7px 10px;
      border-radius: 8px;
      color: #d4deea;
      text-decoration: none;
      font-size: 0.9rem;
    }}
    .nav-link.active {{
      background: linear-gradient(90deg, rgba(77, 157, 247, 0.32), rgba(57, 183, 201, 0.14));
      color: white;
    }}
    .nav-link:hover, .output-link:hover, .action-link:hover {{
      border-color: rgba(77, 157, 247, 0.6);
      background-color: rgba(77, 157, 247, 0.08);
    }}
    .nav-icon {{
      width: 1.35em;
      text-align: center;
      color: var(--accent);
    }}
    .shell {{
      width: min(1480px, calc(100% - 36px));
      margin: 0 auto;
      padding: 18px 0 40px;
    }}
    header {{
      display: grid;
      grid-template-columns: minmax(0, 1fr) auto;
      gap: 20px;
      align-items: end;
      padding: 20px 0 24px;
      border-bottom: 1px solid var(--line);
    }}
    .breadcrumb {{
      color: var(--muted);
      font-size: 0.86rem;
      margin-top: 5px;
    }}
    .brand {{
      color: var(--accent);
      font-size: 0.82rem;
      font-weight: 700;
      letter-spacing: 0.08em;
      text-transform: uppercase;
    }}
    h1 {{
      margin: 6px 0;
      font-size: clamp(1.8rem, 3vw, 3.2rem);
      line-height: 1.05;
      letter-spacing: 0;
    }}
    .subtle {{ color: var(--muted); }}
    .status {{
      display: inline-flex;
      align-items: center;
      gap: 8px;
      padding: 8px 12px;
      border: 1px solid var(--line);
      border-radius: 8px;
      background: rgba(255, 255, 255, 0.04);
      color: var(--muted);
      white-space: nowrap;
    }}
    .dot {{
      width: 9px;
      height: 9px;
      border-radius: 99px;
      background: var(--warn);
    }}
    .dot.ok {{ background: var(--accent-2); }}
    .dot.error {{ background: var(--danger); }}
    .dot.unknown {{ background: var(--warn); }}
    .notice {{
      margin: 20px 0 0;
      padding: 14px 16px;
      border: 1px solid rgba(57, 183, 201, 0.38);
      border-radius: 8px;
      background: rgba(57, 183, 201, 0.09);
      color: #d9f8fc;
    }}
    .notice strong {{
      display: block;
      color: white;
      margin-bottom: 3px;
    }}
    .live-panel {{
      margin: 20px 0 0;
      display: grid;
      grid-template-columns: minmax(0, 1fr) minmax(220px, 0.36fr);
      gap: 14px;
      align-items: stretch;
    }}
    .live-panel .panel-block {{
      margin: 0;
    }}
    .live-status-grid {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(150px, 1fr));
      gap: 10px;
      margin-top: 12px;
    }}
    .live-mini-card {{
      border: 1px solid rgba(148, 163, 184, 0.16);
      border-radius: 8px;
      background: #0b1626;
      padding: 10px;
    }}
    .live-mini-card span {{
      display: block;
      color: var(--muted);
      font-size: 0.75rem;
      margin-bottom: 4px;
    }}
    .live-mini-card strong {{
      overflow-wrap: anywhere;
    }}
    .cards {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
      gap: 14px;
      margin: 22px 0;
    }}
    .card, .plot-card, .output-link, .empty-state, .panel-block, .action-link {{
      border: 1px solid var(--line);
      border-radius: 8px;
      background: linear-gradient(180deg, rgba(19, 35, 56, 0.92), rgba(16, 27, 41, 0.96));
      box-shadow: 0 16px 42px var(--shadow);
    }}
    .card {{
      min-height: 124px;
      padding: 18px;
    }}
    .card .label {{
      color: var(--muted);
      font-size: 0.88rem;
      margin-bottom: 10px;
    }}
    .card .value {{
      font-size: 1.55rem;
      font-weight: 760;
      overflow-wrap: anywhere;
    }}
    .card.good .value {{ color: var(--accent-2); }}
    .card.warn .value {{ color: var(--warn); }}
    .card .detail {{
      color: var(--muted);
      font-size: 0.86rem;
      margin-top: 8px;
      overflow-wrap: anywhere;
    }}
    .status-card {{
      border-color: rgba(57, 183, 201, 0.35);
      background: linear-gradient(180deg, rgba(57, 183, 201, 0.11), rgba(16, 27, 41, 0.96));
    }}
    .section-heading {{
      display: flex;
      justify-content: space-between;
      gap: 16px;
      align-items: baseline;
      margin: 26px 0 12px;
    }}
    h2 {{
      margin: 0;
      font-size: 1.08rem;
      letter-spacing: 0;
    }}
    .plot-grid {{
      display: grid;
      grid-template-columns: repeat(auto-fill, minmax(280px, 1fr));
      grid-auto-rows: 8px;
      grid-auto-flow: dense;
      column-gap: 18px;
      row-gap: 8px;
      align-items: stretch;
    }}
    .plot-card {{
      overflow: hidden;
      min-width: 0;
      min-height: 0;
      padding: 14px;
      display: flex;
      flex-direction: column;
      position: relative;
      grid-column: span var(--col-span, 1);
      grid-row: span var(--row-span, 20);
    }}
    .plot-card.card-sm {{ --col-span: 1; --row-span: 16; }}
    .plot-card.card-md {{ --col-span: 1; --row-span: 20; }}
    .plot-card.card-lg {{ --col-span: 2; --row-span: 38; }}
    .plot-card.card-wide {{ --col-span: 2; --row-span: 24; }}
    .plot-header {{
      display: flex;
      justify-content: space-between;
      gap: 10px;
      align-items: flex-start;
      margin-bottom: 12px;
    }}
    .plot-header h3 {{
      margin: 0;
      font-size: 1rem;
    }}
    .plot-title-group {{
      min-width: 0;
    }}
    .size-controls {{
      display: inline-flex;
      gap: 4px;
      align-items: center;
      flex-wrap: wrap;
      justify-content: flex-end;
    }}
    .size-button, .reset-layout {{
      border: 1px solid rgba(148, 163, 184, 0.28);
      background: rgba(15, 26, 42, 0.88);
      color: var(--muted);
      border-radius: 6px;
      padding: 3px 6px;
      font: inherit;
      font-size: 0.68rem;
      line-height: 1.2;
      cursor: pointer;
    }}
    .size-button:hover, .size-button.active, .reset-layout:hover {{
      color: white;
      border-color: rgba(57, 183, 201, 0.72);
      background: rgba(57, 183, 201, 0.18);
    }}
    .plot-frame {{
      flex: 1 1 auto;
      min-height: 180px;
      display: flex;
      align-items: center;
      justify-content: center;
      overflow: hidden;
      background: #0b1626;
      border-radius: 12px;
      padding: 8px;
      border: 1px solid rgba(148, 163, 184, 0.18);
    }}
    .plot-card.fallback .plot-frame {{
      background: #0b1626;
      padding: 12px;
    }}
    .plot-frame a {{
      display: flex;
      align-items: center;
      justify-content: center;
      width: 100%;
      height: 100%;
      text-decoration: none;
      min-width: 0;
    }}
    .plot-frame img {{
      display: block;
      max-width: 100%;
      max-height: 100%;
      width: 100%;
      height: 100%;
      object-fit: contain;
      border-radius: 6px;
      background: #0f1a2a;
    }}
    .plot-card.fallback .plot-frame img {{
      background: white;
    }}
    .resize-handle {{
      position: absolute;
      right: 8px;
      bottom: 8px;
      width: 16px;
      height: 16px;
      cursor: nwse-resize;
      opacity: 0.75;
      touch-action: none;
    }}
    .resize-handle::before {{
      content: "";
      position: absolute;
      inset: 0;
      background:
        linear-gradient(135deg, transparent 0 45%, rgba(148, 163, 184, 0.9) 46% 52%, transparent 53%),
        linear-gradient(135deg, transparent 0 65%, rgba(148, 163, 184, 0.7) 66% 72%, transparent 73%);
    }}
    .tag {{
      display: inline-flex;
      align-items: center;
      padding: 3px 7px;
      border-radius: 999px;
      border: 1px solid rgba(57, 183, 201, 0.32);
      color: #b9ecf3;
      background: rgba(57, 183, 201, 0.10);
      font-size: 0.72rem;
      white-space: nowrap;
    }}
    .summary-value {{
      color: var(--text);
      font-size: 0.84rem;
      margin-top: 9px;
    }}
    .source {{
      padding: 10px 2px 0;
      color: var(--muted);
      font-size: 0.82rem;
      overflow-wrap: anywhere;
      word-break: break-word;
      line-height: 1.35;
    }}
    .plot-meta {{
      overflow-wrap: anywhere;
      font-size: 0.82rem;
    }}
    .category-label {{
      color: var(--muted);
      font-size: 0.78rem;
      margin-top: 8px;
    }}
    .lower-grid {{
      display: grid;
      grid-template-columns: minmax(280px, 0.95fr) minmax(320px, 1.15fr)
        minmax(280px, 0.95fr) minmax(250px, 0.8fr);
      gap: 14px;
      align-items: start;
      margin-top: 26px;
    }}
    .panel-block {{
      padding: 16px;
    }}
    .panel-block h2 {{
      margin-bottom: 14px;
    }}
    .phase-row {{
      display: grid;
      grid-template-columns: 22px 1fr auto;
      gap: 9px;
      align-items: center;
      padding: 9px 0;
      border-bottom: 1px solid rgba(255, 255, 255, 0.07);
    }}
    .phase-row:last-child {{ border-bottom: 0; }}
    .phase-dot {{
      display: grid;
      place-items: center;
      width: 22px;
      height: 18px;
      border-radius: 99px;
      background: #33465c;
      color: white;
      font-size: 0.58rem;
      font-weight: 800;
    }}
    .phase-row.ok .phase-dot {{ background: var(--accent-2); }}
    .phase-row.error .phase-dot {{ background: var(--danger); }}
    .phase-row.skipped .phase-dot {{ background: var(--warn); }}
    .phase-row.not-run .phase-dot {{ background: #536274; }}
    .phase-detail {{
      color: var(--muted);
      font-size: 0.84rem;
    }}
    table {{
      width: 100%;
      border-collapse: collapse;
      font-size: 0.88rem;
    }}
    th, td {{
      padding: 9px 8px;
      border-bottom: 1px solid rgba(255, 255, 255, 0.08);
      text-align: left;
      vertical-align: top;
    }}
    th {{
      color: var(--muted);
      font-weight: 700;
    }}
    td.num {{ text-align: right; }}
    .table-empty {{
      color: var(--muted);
      text-align: center;
    }}
    .outputs {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(220px, 1fr));
      gap: 10px;
    }}
    .output-link {{
      display: block;
      padding: 9px 10px;
      text-decoration: none;
      box-shadow: none;
    }}
    .output-link strong {{
      display: block;
      margin-bottom: 3px;
    }}
    .output-link span {{
      display: block;
      color: var(--muted);
      font-size: 0.78rem;
      overflow-wrap: anywhere;
      word-break: break-word;
      line-height: 1.35;
    }}
    .outputs-extra {{
      margin-top: 10px;
    }}
    .outputs-extra summary {{
      cursor: pointer;
      color: #b9ecf3;
      font-size: 0.84rem;
      margin-bottom: 10px;
    }}
    .actions {{
      display: grid;
      gap: 10px;
    }}
    .action-link {{
      display: block;
      padding: 12px 13px;
      text-decoration: none;
      box-shadow: none;
    }}
    .action-title {{
      display: block;
      font-weight: 700;
    }}
    .action-subtitle {{
      display: block;
      margin-top: 4px;
      color: var(--muted);
      font-size: 0.9rem;
    }}
    .artifact-path {{
      overflow-wrap: anywhere;
      word-break: break-word;
      line-height: 1.35;
    }}
    .folder-note {{
      margin-top: 12px;
      color: var(--muted);
      font-size: 0.82rem;
      overflow-wrap: anywhere;
    }}
    .empty-state {{
      padding: 20px;
      color: var(--muted);
    }}
    footer {{
      margin-top: 30px;
      color: var(--muted);
      font-size: 0.82rem;
    }}
    @media (max-width: 720px) {{
      .layout {{ grid-template-columns: 1fr; }}
      .sidebar {{
        position: relative;
        height: auto;
        border-right: 0;
        border-bottom: 1px solid var(--line);
      }}
      .shell {{ width: min(100% - 20px, 1480px); padding-top: 14px; }}
      header {{ grid-template-columns: 1fr; align-items: start; }}
      .plot-grid {{ grid-template-columns: 1fr; }}
      .plot-card,
      .plot-card.card-lg,
      .plot-card.card-wide {{ --col-span: 1 !important; }}
      .lower-grid {{ grid-template-columns: 1fr; }}
      .live-panel {{ grid-template-columns: 1fr; }}
    }}
    @media (min-width: 721px) and (max-width: 1320px) {{
      .lower-grid {{ grid-template-columns: repeat(2, minmax(0, 1fr)); }}
    }}
  </style>
</head>
<body>
  <div class="layout">
    <aside class="sidebar">
      {nav_html}
    </aside>
    <main class="shell">
      <header id="dashboard">
        <div>
          <div class="brand">FastMDXplora Results</div>
          <h1>Dashboard</h1>
          <div class="subtle">{escape(title)}</div>
          <div class="breadcrumb">Project / Results / Dashboard</div>
          <div class="subtle">System: {escape(system or "not available")}</div>
        </div>
        <div class="status">
          <span class="dot {escape(status)}"></span>{escape(status.title())} -
          Last generated {escape(generated)}
          <button class="reset-layout" type="button" data-reset-layout>Reset layout</button>
        </div>
      </header>
      {notice_html}
      {live_html}
      <section class="cards" aria-label="Summary metrics">
        {card_html}
      </section>
      {sections_html}
      <section class="lower-grid">
        <div class="panel-block" id="run-status">
          <h2>Run Progress</h2>
          {phase_html}
        </div>
        <div class="panel-block" id="top-metrics">
          <h2>Top Metrics</h2>
          <table>
            <thead>
              <tr><th>Metric</th><th>Average</th><th>Std. Dev.</th><th>Unit</th></tr>
            </thead>
            <tbody>{metrics_html}</tbody>
          </table>
        </div>
        <div class="panel-block" id="recent-outputs">
          <h2>Recent Outputs</h2>
          {link_html}
        </div>
        <div class="panel-block" id="quick-actions">
          <h2>Quick Actions</h2>
          <div class="actions">
            {quick_html}
          </div>
          <div class="folder-note">
            Output folder: {escape(output_folder)}.
            Some browsers block direct local folder opening from static pages.
          </div>
        </div>
      </section>
      <footer>
        Generated by FastMDXplora from recorded run artifacts.
        Missing metrics are omitted.
      </footer>
    </main>
  </div>
  <script>
    (() => {{
      const key = "fastmdx-dashboard-card-layout:" + window.location.pathname;
      const presets = {{
        sm: {{ cols: 1, rows: 16 }},
        md: {{ cols: 1, rows: 20 }},
        lg: {{ cols: 2, rows: 38 }},
        wide: {{ cols: 2, rows: 24 }},
      }};
      let saved = {{}};
      try {{ saved = JSON.parse(localStorage.getItem(key) || "{{}}"); }}
      catch (_) {{ saved = {{}}; }}
      const cards = Array.from(document.querySelectorAll(".plot-card[data-card-key]"));
      const clamp = (value, min, max) => Math.max(min, Math.min(max, value));
      const gridMetrics = (grid) => {{
        const style = getComputedStyle(grid);
        const columns = style.gridTemplateColumns.split(" ").filter(Boolean);
        const firstColumn = parseFloat(columns[0]) || 280;
        return {{
          columns: Math.max(1, columns.length),
          columnWidth: firstColumn,
          columnGap: parseFloat(style.columnGap) || 18,
          rowHeight: parseFloat(style.gridAutoRows) || 8,
          rowGap: parseFloat(style.rowGap) || 8,
        }};
      }};
      const applySpan = (card, cols, rows, persist = true) => {{
        const grid = card.closest(".plot-grid");
        const metrics = grid ? gridMetrics(grid) : {{ columns: 1 }};
        cols = clamp(Math.round(cols || 1), 1, metrics.columns || 1);
        rows = clamp(Math.round(rows || 20), 14, 64);
        card.style.setProperty("--col-span", cols);
        card.style.setProperty("--row-span", rows);
        for (const name of Object.keys(presets)) card.classList.remove("card-" + name);
        for (const button of card.querySelectorAll("[data-card-size]")) {{
          const preset = presets[button.dataset.cardSize];
          button.classList.toggle(
            "active",
            Boolean(preset && preset.cols === cols && preset.rows === rows)
          );
        }}
        if (persist) {{
          saved[card.dataset.cardKey] = {{ cols, rows }};
          try {{ localStorage.setItem(key, JSON.stringify(saved)); }}
          catch (_) {{}}
        }}
      }};
      const applySize = (card, size, persist = true) => {{
        const preset = presets[size] || presets.md;
        applySpan(card, preset.cols, preset.rows, persist);
      }};
      for (const card of cards) {{
        const restored = saved[card.dataset.cardKey];
        if (restored && Number.isFinite(restored.cols) && Number.isFinite(restored.rows)) {{
          applySpan(card, restored.cols, restored.rows, false);
        }} else {{
          applySize(card, "md", false);
        }}
        for (const button of card.querySelectorAll("[data-card-size]")) {{
          button.addEventListener("click", () => applySize(card, button.dataset.cardSize));
        }}
        const handle = card.querySelector(".resize-handle");
        if (handle) {{
          handle.addEventListener("pointerdown", (event) => {{
            event.preventDefault();
            handle.setPointerCapture(event.pointerId);
            const grid = card.closest(".plot-grid");
            const metrics = gridMetrics(grid);
            const start = {{
              x: event.clientX,
              y: event.clientY,
              cols: parseInt(getComputedStyle(card).getPropertyValue("--col-span"), 10) || 1,
              rows: parseInt(getComputedStyle(card).getPropertyValue("--row-span"), 10) || 20,
            }};
            const move = (moveEvent) => {{
              const colUnit = metrics.columnWidth + metrics.columnGap;
              const rowUnit = metrics.rowHeight + metrics.rowGap;
              const cols = start.cols + Math.round((moveEvent.clientX - start.x) / colUnit);
              const rows = start.rows + Math.round((moveEvent.clientY - start.y) / rowUnit);
              applySpan(card, cols, rows, false);
            }};
            const done = () => {{
              handle.removeEventListener("pointermove", move);
              handle.removeEventListener("pointerup", done);
              handle.removeEventListener("pointercancel", done);
              const cols = parseInt(getComputedStyle(card).getPropertyValue("--col-span"), 10) || 1;
              const rows = parseInt(getComputedStyle(card).getPropertyValue("--row-span"), 10) || 20;
              applySpan(card, cols, rows, true);
            }};
            handle.addEventListener("pointermove", move);
            handle.addEventListener("pointerup", done);
            handle.addEventListener("pointercancel", done);
          }});
        }}
      }}
      const reset = document.querySelector("[data-reset-layout]");
      if (reset) {{
        reset.addEventListener("click", () => {{
          saved = {{}};
          try {{ localStorage.removeItem(key); }} catch (_) {{}}
          for (const card of cards) applySize(card, "md", false);
        }});
      }}
    }})();
  </script>
</body>
</html>
"""


def _render_card(card: DashboardCard) -> str:
    detail = f'<div class="detail">{escape(card.detail)}</div>' if card.detail else ""
    classes = f"card {card.kind}"
    if card.label == "Project status":
        classes += " status-card"
    return (
        f'<article class="{escape(classes)}">'
        f'<div class="label">{escape(card.label)}</div>'
        f'<div class="value">{escape(card.value)}</div>'
        f"{detail}"
        "</article>"
    )


def _render_section(section: DashboardSection) -> str:
    panels = "\n".join(_render_panel(panel) for panel in section.panels)
    count = f"{len(section.panels)} artifact" + ("" if len(section.panels) == 1 else "s")
    classes = f"analysis-section section-{section.anchor}"
    return (
        f'<section class="{escape(classes)}" id="{escape(section.anchor)}">'
        '<div class="section-heading">'
        f"<h2>{escape(section.title)}</h2>"
        f'<span class="subtle">{escape(count)}</span>'
        "</div>"
        '<div class="plot-grid">'
        f"{panels}"
        "</div>"
        "</section>"
    )


def _render_panel(panel: DashboardPanel) -> str:
    card_class = (
        "plot-card card-md"
        if panel.mode == "dashboard view"
        else "plot-card fallback card-md"
    )
    summary = (
        f'<div class="summary-value">{escape(panel.summary)}</div>'
        if panel.summary else ""
    )
    category = (
        f'<div class="category-label">{escape(panel.category)}</div>'
        if panel.category else ""
    )
    original = ""
    if panel.original_source != panel.source:
        original = (
            f' - original: <a href="{escape(panel.original_href)}">'
            f"{escape(panel.original_source)}</a>"
        )
    card_key = _anchor(panel.original_source)
    return (
        f'<article class="{escape(card_class)}" id="{escape(_anchor(panel.title))}" '
        f'data-card-key="{escape(card_key)}">'
        '<div class="plot-header">'
        '<div class="plot-title-group">'
        f"<h3>{escape(panel.title)}</h3>"
        f'<span class="tag">{escape(panel.mode)}</span>'
        "</div>"
        '<div class="size-controls" aria-label="Card size">'
        '<button type="button" class="size-button" data-card-size="sm">S</button>'
        '<button type="button" class="size-button active" data-card-size="md">M</button>'
        '<button type="button" class="size-button" data-card-size="lg">L</button>'
        '<button type="button" class="size-button" data-card-size="wide">Wide</button>'
        "</div>"
        "</div>"
        '<div class="plot-frame">'
        f'<a href="{escape(panel.href)}">'
        f'<img src="{escape(panel.href)}" alt="{escape(panel.title)} plot">'
        "</a>"
        "</div>"
        f"{summary}"
        f"{category}"
        f'<div class="source plot-meta artifact-path">{escape(panel.source)}{original}</div>'
        '<div class="resize-handle" title="Drag to resize"></div>'
        "</article>"
    )


def _render_link(link: DashboardLink) -> str:
    detail = escape(link.detail) if link.detail else escape(link.href)
    return (
        f'<a class="output-link" href="{escape(link.href)}">'
        f'<strong class="output-title">{escape(link.label)}</strong>'
        f'<span class="output-subtitle artifact-path">{detail}</span>'
        "</a>"
    )


def _render_output_list(links: list[DashboardLink]) -> str:
    visible = links[:10]
    hidden = links[10:]
    visible_html = "\n".join(_render_link(link) for link in visible)
    html = f'<div class="outputs output-list">{visible_html}</div>'
    if hidden:
        hidden_html = "\n".join(_render_link(link) for link in hidden)
        html += (
            '<details class="outputs-extra">'
            f"<summary>Show all outputs ({len(hidden)} more)</summary>"
            f'<div class="outputs output-list">{hidden_html}</div>'
            "</details>"
        )
    return html


def _render_static_live_panel(project_root: Path) -> str:
    from fastmdxplora.live.telemetry import analyze_health, read_metrics, read_status

    serve_command = f"fastmdx dashboard serve --output {project_root.as_posix()}"
    escaped_command = escape(serve_command)
    status = read_status(project_root)
    metrics = read_metrics(project_root)
    health = analyze_health(status, metrics)
    if not status:
        body = (
            "Live simulation telemetry was not recorded for this run. Start "
            f"<code>{escaped_command}</code> during a "
            "simulation to monitor progress."
        )
        stage = "not available"
        updated = "not available"
        platform = "not available"
    else:
        body = escape(str(health.get("explanation") or "not available"))
        stage = escape(str(status.get("stage") or "not available"))
        updated = escape(str(status.get("last_update_timestamp") or "not available"))
        platform = escape(str(status.get("platform") or "not available"))
    health_state = escape(str(health.get("state") or "unknown"))
    health_message = escape(str(health.get("message") or "Live telemetry is not available."))
    return (
        '<section class="live-panel" id="live-simulation">'
        '<div class="panel-block">'
        "<h2>Live Simulation</h2>"
        f"<p>{body}</p>"
        '<div class="live-status-grid">'
        f'<div class="live-mini-card"><span>Status</span><strong>{health_state}</strong></div>'
        f'<div class="live-mini-card"><span>Stage</span><strong>{stage}</strong></div>'
        f'<div class="live-mini-card"><span>Platform</span><strong>{platform}</strong></div>'
        f'<div class="live-mini-card"><span>Last update</span><strong>{updated}</strong></div>'
        "</div>"
        "</div>"
        '<div class="panel-block">'
        "<h2>Monitoring</h2>"
        f"<p>{health_message}</p>"
        '<p class="subtle">For real-time charts and events, run '
        f"<code>{escaped_command}</code>.</p>"
        "</div>"
        "</section>"
    )


def _render_sidebar(sections: list[DashboardSection], links: list[DashboardLink]) -> str:
    analysis_links = {section.title: f"#{section.anchor}" for section in sections}
    report_links = {link.label: link.href for link in links}
    parts = [
        '<div class="logo">',
        '<div class="mark">FX</div>',
        '<div><div class="logo-title">FastMDXplora</div>',
        '<div class="logo-subtitle">Explore. Analyze. Visualize. Share.</div></div>',
        '</div>',
        '<nav aria-label="Dashboard navigation">',
        '<div class="nav-section"><p class="nav-heading">Overview</p>',
        _nav_link("#dashboard", "Dashboard", "active", "D"),
        _nav_link("#live-simulation", "Live Simulation", "", "L"),
        _nav_link("#top-metrics", "System Info", "", "I"),
        _nav_link("#run-status", "Run Status", "", "S"),
        '</div>',
        '<div class="nav-section"><p class="nav-heading">Analysis</p>',
    ]
    for label in SECTION_ORDER:
        href = analysis_links.get(label)
        if href:
            parts.append(_nav_link(href, label, "", "-"))
    parts.extend(
        [
            '</div>',
            '<div class="nav-section"><p class="nav-heading">Reports</p>',
        ]
    )
    for label in ("Markdown report", "Slide deck", "Dashboard HTML", "Project bundle"):
        href = report_links.get(label)
        if href:
            display = {
                "Markdown report": "Markdown Report",
                "Slide deck": "Slides PPTX",
                "Dashboard HTML": "Dashboard HTML",
                "Project bundle": "Bundle",
            }[label]
            parts.append(_nav_link(href, display, "", "R"))
    parts.extend(["</div>", "</nav>"])
    return "\n".join(parts)


def _nav_link(href: str, label: str, classes: str, icon: str) -> str:
    class_attr = f"nav-link {classes}".strip()
    return (
        f'<a class="{escape(class_attr)}" href="{escape(href)}">'
        f'<span class="nav-icon">{escape(icon)}</span>{escape(label)}</a>'
    )


def _render_phase_row(row: PhaseRow) -> str:
    symbol = {
        "ok": "OK",
        "error": "!",
        "skipped": "-",
        "not-run": "-",
    }.get(row.status, "?")
    return (
        f'<div class="phase-row {escape(row.status)}">'
        f'<span class="phase-dot">{escape(symbol)}</span>'
        f"<span>{escape(row.name)}</span>"
        f'<span class="phase-detail">{escape(row.detail)}</span>'
        "</div>"
    )


def _render_metric_row(row: MetricRow) -> str:
    return (
        "<tr>"
        f"<td>{escape(row.metric)}</td>"
        f'<td class="num">{escape(row.average)}</td>'
        f'<td class="num">{escape(row.stddev)}</td>'
        f"<td>{escape(row.unit)}</td>"
        "</tr>"
    )


def _quick_action_links(links: list[DashboardLink]) -> list[DashboardLink]:
    by_label = {link.label: link for link in links}
    actions: list[DashboardLink] = []
    for label, display in (
        ("Markdown report", "Open Markdown Report"),
        ("Slide deck", "Open Slides"),
        ("Project bundle", "Open Bundle"),
        ("Analysis manifest", "Open Analysis Manifest"),
        ("Dashboard HTML", "Open Dashboard HTML"),
    ):
        link = by_label.get(label)
        if link:
            actions.append(DashboardLink(display, link.href, link.detail))
    return actions


def _render_quick_action(link: DashboardLink) -> str:
    return (
        f'<a class="action-link" href="{escape(link.href)}">'
        f'<span class="action-title">{escape(link.label)}</span>'
        f'<span class="action-subtitle">{escape(link.detail)}</span>'
        "</a>"
    )


def _anchor(value: str) -> str:
    chars = []
    for char in value.lower():
        if char.isalnum():
            chars.append(char)
        elif chars and chars[-1] != "-":
            chars.append("-")
    anchor = "".join(chars).strip("-")
    return anchor or "section"
