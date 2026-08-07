"""Helpers for making report artifacts reflect the phases that actually ran."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class PhaseContext:
    """Phase presence/status loaded from the project manifest."""

    setup_present: bool = False
    simulation_present: bool = False
    analysis_present: bool = False
    report_present: bool = False

    @property
    def is_full_pipeline(self) -> bool:
        return self.setup_present and self.simulation_present and self.analysis_present

    @property
    def is_analysis_from_existing_trajectory(self) -> bool:
        return self.analysis_present and not self.setup_present and not self.simulation_present


def load_phase_context(project_root: Path) -> PhaseContext:
    """Return the phases recorded in ``manifest.json``.

    Older or partial outputs may not have a manifest; in that case, fall back
    to the presence of phase artifact directories so report-only mode still
    produces useful wording.
    """
    phases: set[str] = set()
    manifest_path = project_root / "manifest.json"
    try:
        with manifest_path.open(encoding="utf-8") as fh:
            manifest = json.load(fh)
        raw_phases = manifest.get("phases", [])
        if isinstance(raw_phases, list):
            phases.update(
                str(phase.get("name", ""))
                for phase in raw_phases
                if isinstance(phase, dict) and phase.get("status") != "skipped"
            )
    except (FileNotFoundError, OSError, json.JSONDecodeError):
        phases = set()

    if not phases:
        if (project_root / "setup").is_dir():
            phases.add("setup")
        if (project_root / "simulation").is_dir():
            phases.add("simulation")
        if (project_root / "analysis").is_dir():
            phases.add("analysis")
        if (project_root / "report").is_dir():
            phases.add("report")

    return PhaseContext(
        setup_present="setup" in phases,
        simulation_present="simulation" in phases,
        analysis_present="analysis" in phases,
        report_present="report" in phases,
    )


def _system_label(system: object) -> str:
    """What to call the system in a title.

    The name of the structure, not the path to it. A report is a document
    somebody sends to a colleague or attaches to a submission, and the title
    read "FastMDXplora Study --- /home/claude/ala3.pdb": the author's home
    directory, on the first line, in a file meant to leave the machine. A
    four-character PDB entry is already a name and is kept as one.
    """
    from pathlib import Path

    text = str(system)
    if len(text) == 4 and text.isalnum():
        return text.upper()
    stem = Path(text).name
    # Repeatedly, because a compressed structure carries two: taking one left
    # "trp_cage.pdb" as the name of a study.
    changed = True
    while changed:
        changed = False
        for suffix in (".gz", ".bz2", ".pdb", ".cif", ".pdbx", ".ent"):
            if stem.lower().endswith(suffix) and len(stem) > len(suffix):
                stem = stem[: -len(suffix)]
                changed = True
    return stem or text
