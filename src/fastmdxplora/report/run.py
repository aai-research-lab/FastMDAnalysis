"""Top-level report phase entry point.

Combines outputs from the document, slides, and bundle generators into a
unified report artifact set.
"""

from __future__ import annotations

from fastmdxplora.utils.logging import get_logger
from pathlib import Path
from typing import TYPE_CHECKING, Any

from fastmdxplora.report.bundle import build_bundle
from fastmdxplora.gui.report_dashboard import build_dashboard
from fastmdxplora.report.document import build_document
from fastmdxplora.report.region_highlights import build_region_highlight_artifacts
from fastmdxplora.report.slides import build_slides
from fastmdxplora.report.summary_figure import build_analysis_summary_figure

if TYPE_CHECKING:
    from fastmdxplora.orchestrator import FastMDXplora

logger = get_logger("report")


DEFAULTS: dict[str, Any] = {
    "document": True,
    "slides": True,
    "bundle": True,
    "title": None,  # auto-derive from system
    "author": None,
    "include_methods": True,
    "include_reproducibility": True,
    "region_highlights": None,
}


def run(
    *,
    orchestrator: "FastMDXplora",
    output_dir: Path,
    **options: Any,
) -> list[str]:
    """Run the report phase.

    Parameters
    ----------
    orchestrator : FastMDXplora
        Parent orchestrator (used to locate sibling phase outputs).
    output_dir : Path
        Where to write report artifacts.
    **options
        See ``DEFAULTS``.

    Returns
    -------
    list of str
        Paths (relative to ``output_dir``) of artifacts produced.
    """
    params = {**DEFAULTS, **options}
    title = params["title"] or f"FastMDXplora Study — {orchestrator.system}"

    presenter = getattr(orchestrator, "_presenter", None)
    artifacts: list[str] = []

    region_artifacts = build_region_highlight_artifacts(
        project_root=orchestrator.output_dir,
        output_dir=output_dir,
        region_highlights=params.get("region_highlights"),
    )
    artifacts.extend(region_artifacts)
    if presenter:
        for art in region_artifacts:
            presenter.step(f"Wrote {art}")

    summary_artifacts = build_analysis_summary_figure(
        project_root=orchestrator.output_dir,
        output_dir=output_dir,
    )
    artifacts.extend(summary_artifacts)
    if presenter:
        for art in summary_artifacts:
            presenter.step(f"Wrote {art}")

    if params["document"]:
        doc_artifacts = build_document(
            orchestrator=orchestrator,
            output_dir=output_dir,
            title=title,
            author=params["author"],
            include_methods=params["include_methods"],
            include_reproducibility=params["include_reproducibility"],
        )
        artifacts.extend(doc_artifacts)
        if presenter:
            presenter.step("Wrote report.md")

    if params["slides"]:
        slide_artifacts = build_slides(
            orchestrator=orchestrator,
            output_dir=output_dir,
            title=title,
        )
        artifacts.extend(slide_artifacts)
        if presenter:
            for art in slide_artifacts:
                presenter.step(f"Wrote {art}")

    dashboard_artifacts = build_dashboard(
        orchestrator=orchestrator,
        output_dir=output_dir,
        title=title,
        include_bundle_link=bool(params["bundle"]),
    )
    artifacts.extend(dashboard_artifacts)
    if presenter:
        for art in dashboard_artifacts:
            presenter.step(f"Wrote {art}")

    if params["bundle"]:
        bundle_artifacts = build_bundle(
            orchestrator=orchestrator,
            output_dir=output_dir,
        )
        artifacts.extend(bundle_artifacts)
        if presenter:
            presenter.step("Wrote project_bundle.zip")

    logger.debug("report: wrote %d artifact(s) to %s", len(artifacts), output_dir)
    return artifacts
