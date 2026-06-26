"""Structured study report (Markdown).

Produces a publication-style report from the project state:

  - Header (title, authors, date)
  - Methods (auto-populated from setup + simulation parameter manifests)
  - Results (figures + summary tables from the analysis manifest)
  - Discussion (stub for the user to fill in)
  - Citation (FastMDXplora + JCC paper)
  - Reproducibility appendix (command-line invocation, software versions,
    parameter manifests, input hashes)

PDF rendering of this report is an optional add-on (requires extra
dependencies); the Markdown source is always produced.
"""

from __future__ import annotations

import json
from fastmdxplora.utils.logging import get_logger
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import TYPE_CHECKING
from urllib.parse import quote

from fastmdxplora.report.context import PhaseContext, load_phase_context
from fastmdxplora.analysis.energy import ensure_energy_report_assets
from fastmdxplora.report.captions import caption_for_figure

if TYPE_CHECKING:
    from fastmdxplora.orchestrator import FastMDXplora

logger = get_logger("report.document")


def _one_line(value: object, *, limit: int = 1000) -> str:
    text = str(value)
    text = " ".join(text.replace("\t", " ").splitlines())
    text = " ".join(text.split())
    if len(text) > limit:
        return text[: limit - 1].rstrip() + "..."
    return text


def _md_text(value: object, *, limit: int = 1000) -> str:
    text = _one_line(value, limit=limit)
    for char in "\\`*_{}[]()#+-.!|<>":
        text = text.replace(char, f"\\{char}")
    return text


def _code_text(value: object, *, limit: int = 1000) -> str:
    return _one_line(value, limit=limit).replace("`", "'")


def _link_target(path: str) -> str:
    return quote(path, safe="/._-")


def _load_json_safely(path: Path) -> dict | None:
    try:
        with path.open(encoding="utf-8") as fh:
            return json.load(fh)
    except FileNotFoundError:
        return None
    except json.JSONDecodeError:
        logger.warning("Could not parse JSON manifest at %s", path)
        return None




def _simulation_energy_section(project_root: Path) -> list[str]:
    """Generate and embed simulation energy diagnostics when energy.csv exists."""
    energy_csv = project_root / "simulation" / "energy.csv"
    if not energy_csv.exists():
        return []

    lines: list[str] = ["### Simulation Energy Diagnostics", ""]
    lines.append(
        "FastMDXplora generated diagnostic plots from "
        "`simulation/energy.csv` to summarize simulation stability."
    )
    lines.append("")

    try:
        assets = ensure_energy_report_assets(energy_csv, project_root / "simulation")
    except Exception as exc:  # noqa: BLE001
        logger.warning("Could not generate energy diagnostics: %s", exc)
        lines.append(
            "_Energy CSV was found, but diagnostic plots could not be generated._"
        )
        lines.append("")
        return lines

    for key in ("energy_trace", "simulation_health"):
        fig = assets.get(key)
        if fig is None or not fig.exists():
            continue

        rel = fig.relative_to(project_root).as_posix()
        alt_text = _md_text(key.replace("_", " ").title())
        caption = _md_text(caption_for_figure("simulation", fig.stem), limit=2000)

        lines.append(f"![{alt_text}]({_link_target(rel)})")
        lines.append("")
        lines.append(f"**Figure:** {caption}")
        lines.append("")

    summary = assets.get("energy_summary")
    if summary is not None and summary.exists():
        rel = summary.relative_to(project_root).as_posix()
        lines.append(f"Energy summary table: [`simulation/energy_summary.csv`]({_link_target(rel)})")
        lines.append("")

    return lines

def _summary_section(phase_context: PhaseContext) -> str:
    if phase_context.is_full_pipeline:
        summary = (
            "This report was generated automatically by FastMDXplora from the "
            "outputs of an end-to-end molecular dynamics study."
        )
    elif phase_context.is_analysis_from_existing_trajectory:
        summary = (
            "This report was generated from an existing trajectory. Setup and "
            "simulation were not run in this workflow."
        )
    elif phase_context.analysis_present:
        summary = (
            "This report summarizes the FastMDXplora phases recorded for this "
            "workflow."
        )
    else:
        summary = (
            "This report summarizes the available FastMDXplora outputs for "
            "this workflow."
        )
    return f"## Summary\n\n{summary}"


def _methods_section(project_root: Path, phase_context: PhaseContext) -> str:
    setup = _load_json_safely(project_root / "setup" / "setup_parameters.json") or {}
    sim = _load_json_safely(project_root / "simulation" / "simulation_parameters.json") or {}
    setup_params = setup.get("parameters", {})
    sim_params = sim.get("parameters", {})

    lines = ["## Methods", ""]
    lines.append("### System preparation")
    if setup_params:
        lines.append("")
        lines.append(
            "The input system was prepared using FastMDXplora's automated "
            "setup pipeline with the following parameters:"
        )
        lines.append("")
        for k, v in setup_params.items():
            lines.append(f"- **{_md_text(k)}**: `{_code_text(v)}`")
    else:
        lines.append("")
        if phase_context.setup_present:
            lines.append("Setup ran in this workflow, but parameters were not recorded.")
        else:
            lines.append("Setup was not run in this workflow.")

    lines.append("")
    lines.append("### Molecular dynamics simulation")
    if sim_params:
        lines.append("")
        lines.append(
            "Production MD was performed with the following simulation parameters:"
        )
        lines.append("")
        for k, v in sim_params.items():
            lines.append(f"- **{_md_text(k)}**: `{_code_text(v)}`")
    else:
        lines.append("")
        if phase_context.simulation_present:
            lines.append(
                "Simulation ran in this workflow, but parameters were not recorded."
            )
        elif phase_context.analysis_present:
            lines.append(
                "Simulation was not run in this workflow. Analysis was performed "
                "on externally provided or previously generated trajectory/topology "
                "files."
            )
        else:
            lines.append("Simulation was not run in this workflow.")

    return "\n".join(lines)


def _results_section(project_root: Path) -> str:
    analysis_manifest = _load_json_safely(
        project_root / "analysis" / "analysis_manifest.json"
    ) or {}
    plan: list[str] = analysis_manifest.get("plan", [])
    results = analysis_manifest.get("results", {})

    lines = ["## Results", ""]
    energy_lines = _simulation_energy_section(project_root)

    if not plan:
        lines.append("No analyses were executed in this session.")
        if energy_lines:
            lines.append("")
            lines.extend(energy_lines)
        return "\n".join(lines)

    n_frames = analysis_manifest.get("n_frames")
    n_residues = analysis_manifest.get("n_residues")
    if n_frames is not None and n_residues is not None:
        lines.append(
            f"Analysis was performed on a trajectory of {n_frames} frames "
            f"and {n_residues} residues."
        )
    lines.append(f"Analyses performed: {', '.join(_md_text(a) for a in plan)}.")
    lines.append("")

    for analysis in plan:
        # Pretty heading: uppercase short names, title-case longer ones
        heading = analysis.upper() if len(analysis) <= 4 else analysis.title()
        heading = _md_text(heading)
        lines.append(f"### {heading}")
        lines.append("")

        # Per-analysis result row from the analysis manifest
        result_meta = results.get(analysis, {})
        status = result_meta.get("status", "unknown")
        if status != "ok":
            lines.append(
                f"_This analysis did not complete successfully (status: "
                f"`{status}`)._"
            )
            if result_meta.get("message"):
                lines.append(f"Reason: {_md_text(result_meta['message'])}")
            lines.append("")
            continue

        # Options come from each analysis's own options.json (more reliable
        # than the manifest because the per-analysis file records the
        # actual fully-resolved options after defaults are applied).
        opts_file = project_root / "analysis" / analysis / "options.json"
        per_analysis = _load_json_safely(opts_file) or {}
        opts = per_analysis.get("options", {})
        selection = per_analysis.get("selection")

        if opts or selection:
            lines.append("**Parameters:**")
            if selection:
                lines.append(f"- `selection`: `{_code_text(selection)}`")
            for k, v in opts.items():
                lines.append(f"- `{_code_text(k)}`: `{_code_text(v)}`")
        else:
            lines.append("_Ran with default options._")
        lines.append("")

        # Embed all figures in the analysis directory. Multi-method
        # analyses (cluster, dimred) emit several PNGs; sort them so the
        # report renders deterministically.
        figs_dir = project_root / "analysis" / analysis
        figures = sorted(figs_dir.glob("*.png")) if figs_dir.exists() else []
        if figures:
            for fig in figures:
                # Markdown/HTML image links require forward slashes on every
                # OS; str(WindowsPath) would emit backslashes and break them.
                rel = fig.relative_to(project_root).as_posix()
                alt_text = _md_text(f"{analysis} — {fig.stem}")
                caption = _md_text(caption_for_figure(analysis, fig.stem), limit=2000)
                lines.append(f"![{alt_text}]({_link_target(rel)})")
                lines.append("")
                lines.append(f"**Figure:** {caption}")
                lines.append("")
        else:
            lines.append("_No figure was produced for this analysis._")
            lines.append("")

    if energy_lines:
        lines.append("")
        lines.extend(energy_lines)

    return "\n".join(lines)


def _citation_section() -> str:
    from fastmdxplora import __citation__

    return "\n".join(
        [
            "## Citation",
            "",
            "If you use FastMDXplora in your work, please cite:",
            "",
            f"> {__citation__}",
            "",
            "BibTeX:",
            "",
            "```bibtex",
            "@article{aina2026fastmd,",
            "  author  = {Aina, Adekunle and Kwan, Derrick},",
            "  title   = {FastMDAnalysis: Software for Automated Analysis of "
            "Molecular Dynamics Trajectories},",
            "  journal = {Journal of Computational Chemistry},",
            "  volume  = {47},",
            "  number  = {8},",
            "  pages   = {e70350},",
            "  year    = {2026},",
            "  doi     = {10.1002/jcc.70350},",
            "}",
            "```",
        ]
    )


def _reproducibility_section(
    orchestrator: "FastMDXplora",
    phase_context: PhaseContext,
) -> str:
    from fastmdxplora import __version__

    lines = ["## Reproducibility", ""]
    lines.append(f"- **FastMDXplora version**: `{__version__}`")
    lines.append(f"- **Python**: `{sys.version.split()[0]}`")
    lines.append(f"- **Platform**: `{platform.platform()}`")
    lines.append(f"- **System input**: `{_code_text(orchestrator.system)}`")
    lines.append(f"- **Output directory**: `{_code_text(orchestrator.output_dir)}`")
    lines.append("")
    manifests: list[str] = []
    if phase_context.setup_present:
        manifests.append("`setup/setup_parameters.json`")
    if phase_context.simulation_present:
        manifests.append("`simulation/simulation_parameters.json`")
    if phase_context.analysis_present:
        manifests.append("`analysis/analysis_manifest.json`")
    if manifests:
        lines.append(
            "Per-phase parameter manifests for phases in this workflow are "
            f"preserved at {', '.join(manifests)}. The complete session manifest "
            "is at `manifest.json` at the project root."
        )
    else:
        lines.append(
            "The complete session manifest is at `manifest.json` at the project root."
        )
    return "\n".join(lines)


def build_document(
    *,
    orchestrator: "FastMDXplora",
    output_dir: Path,
    title: str,
    author: str | None,
    include_methods: bool,
    include_reproducibility: bool,
) -> list[str]:
    """Render the Markdown study report.

    Returns
    -------
    list of str
        Artifact paths relative to ``output_dir``.
    """
    project_root = orchestrator.output_dir
    phase_context = load_phase_context(project_root)
    now = datetime.now(timezone.utc).strftime("%Y-%m-%d")

    sections: list[str] = []
    header = [f"# {_md_text(title, limit=200)}", ""]
    if author:
        header.append(f"_Author: {_md_text(author, limit=200)}_  ")
    header.append(f"_Generated: {now} (UTC)_  ")
    header.append(f"_Tool: FastMDXplora_")
    sections.append("\n".join(header))

    sections.append(_summary_section(phase_context))

    if include_methods:
        sections.append(_methods_section(project_root, phase_context))

    sections.append(_results_section(project_root))

    sections.append(
        "## Discussion\n\n_This section is intended for the user to complete. "
        "FastMDXplora provides the analytical scaffolding; scientific "
        "interpretation remains the researcher's responsibility._"
    )

    sections.append(_citation_section())

    if include_reproducibility:
        sections.append(_reproducibility_section(orchestrator, phase_context))

    doc = "\n\n".join(sections) + "\n"
    md_path = output_dir / "report.md"
    md_path.write_text(doc, encoding="utf-8")

    return ["report.md"]
