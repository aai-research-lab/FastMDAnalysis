"""Retained stage coordinates, read without starting a simulation engine."""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any
from xml.etree import ElementTree

from fastmdxplora.record_access import confined_records, is_file, open_text, stat
from fastmdxplora.refusals import OutsideWorkspace


MAX_SOURCE_BYTES = 32 * 1024 * 1024
STAGES = (
    ("input", "Structure", ("setup/input.pdb",), None),
    ("prepared", "Prepared", ("setup/prepared.pdb",), None),
    ("solvated", "Solvated", ("setup/solvated.pdb", "setup/topology.pdb"), None),
    ("minimized", "Minimized", ("simulation/state_minimized.xml",), "simulation/topology.pdb"),
    ("equilibrated", "Equilibrated", (), None),
    ("production", "Production", ("simulation/state_final.xml",), "simulation/topology.pdb"),
)


def _public(root: Path, path: Path) -> bool:
    from fastmdxplora.gui.server import _private

    try:
        return not _private(path.relative_to(root).parts)
    except ValueError:
        return False


def _text(path: Path, limit: int = MAX_SOURCE_BYTES) -> str | None:
    if stat(path).st_size > limit:
        return None
    with open_text(path, encoding="utf-8") as handle:
        text = handle.read(limit + 1)
    return text if len(text) <= limit else None


def _record(root: Path, relative: str) -> dict[str, Any]:
    try:
        value = json.loads(_text(root / relative, 1_000_000) or "{}")
        return value if isinstance(value, dict) else {}
    except (OSError, ValueError, UnicodeError):
        return {}


def _version(root: Path, sources: list[str]) -> str:
    parts = []
    for source in sources:
        info = stat(root / source)
        parts.append(f"{source}:{info.st_mtime_ns}:{info.st_size}")
    return hashlib.sha256("|".join(parts).encode()).hexdigest()[:20]


def _state(value: Any) -> str:
    word = str(value or "").lower()
    if word in {"completed", "complete", "ok"}:
        return "completed"
    if word in {"current", "running"}:
        return "current"
    if word in {"failed", "error", "interrupted", "stopped", "skipped"}:
        return word
    return "unknown"


def _entries(root: Path) -> list[dict[str, Any]]:
    from fastmdxplora.gui.telemetry import status_as_it_stands

    status = status_as_it_stands(root, include_campaign=False)
    states = status.get("stage_states")
    states = states if isinstance(states, dict) else {}
    manifest = _record(root, "manifest.json")
    phases = {str(p.get("name") or p.get("phase") or "") for p in manifest.get("phases", [])
              if isinstance(p, dict)} if isinstance(manifest.get("phases"), list) else set()
    has_setup = "setup" in phases or any(is_file(root / p) for p in
        ("setup/input.pdb", "setup/prepared.pdb", "setup/solvated.pdb", "setup/topology.pdb"))
    has_simulation = "simulation" in phases or str(status.get("stage", "")).lower() in {
        "minimization", "nvt", "npt", "production"} or any(is_file(root / p) for p in
        ("simulation/state_minimized.xml", "simulation/state_final.xml"))
    answer = []
    for key, label, candidates, topology in STAGES:
        if not (has_setup if key in {"input", "prepared", "solvated"} else has_simulation):
            continue
        source = next((p for p in candidates if is_file(root / p)), None)
        sources = [source] if source else []
        if topology:
            if not is_file(root / topology):
                source = None
            else:
                sources.append(topology)
        state = _state(states.get({"minimized": "minimization"}.get(key, key)))
        if key == "equilibrated":
            equi = [_state(states.get(name)) for name in ("nvt", "npt")]
            if "current" in equi:
                state = "current"
            elif all(s in {"completed", "skipped"} for s in equi):
                state = "completed" if "completed" in equi else "skipped"
            else:
                state = next((s for s in equi if s in {"failed", "error", "stopped", "interrupted"}), "unknown")
        if key in {"input", "prepared", "solvated"} and source:
            state = "retained"
        current = str(status.get("stage", "")).lower()
        matches = current == {"minimized": "minimization"}.get(key, key)
        matches = matches or (key == "equilibrated" and current in {"nvt", "npt"})
        if matches and _state(status.get("status")) in {"completed", "current", "failed", "error", "stopped", "interrupted"}:
            state = _state(status["status"])
        available = bool(source)
        reason = ""
        if available and any(stat(root / p).st_size > MAX_SOURCE_BYTES for p in sources):
            available = False
            reason = "The retained coordinates exceed the display size limit."
        if not available and not reason:
            reason = ("A separate end-of-equilibration coordinate file was not retained."
                      if key == "equilibrated" else "No retained coordinates and matching topology are available for this stage.")
        entry = {"id": key, "label": label, "state": state, "available": available,
                 "source": source or "", "topology": topology if source and topology else "",
                 "reason": reason, "version": _version(root, sources) if available else ""}
        answer.append(entry)
    return answer


def timeline_payload(root: str | Path) -> dict[str, Any]:
    """Name retained files separately from telemetry's completion claims."""
    root = Path(root).resolve()
    try:
        with confined_records(root, lambda path: _public(root, path)) as scope:
            entries = _entries(root)
        if scope.denied:
            return {"available": False, "stages": [], "reason": "A stage source could not be confined to the study's public files."}
        return {"available": True, "stages": entries, "reason": ""}
    except (OSError, OutsideWorkspace):
        return {"available": False, "stages": [], "reason": "Stage records could not be read within the study's public files."}


def state_coordinates(topology: str, xml: str) -> tuple[str | None, str]:
    """OpenMM State version 1 stores positions and cell vectors in nm.

    The topology copied beside that State supplies atom order and identity.
    PDB coordinates and cell edges are in angstroms. No Context or System
    is created, and velocities, energies and forces are never interpreted.
    """
    if "<!DOCTYPE" in xml.upper() or "<!ENTITY" in xml.upper():
        return None, "XML declarations and entities are not accepted in a retained State."
    try:
        state = ElementTree.fromstring(xml)
        if state.tag != "State" or state.get("type") != "State" or state.get("version") != "1":
            return None, "This is not a supported OpenMM State version."
        positions = state.find("Positions")
        if positions is None or any(p.tag != "Position" for p in positions):
            return None, "The retained State has no readable positions."
        xyz = [[float(p.attrib[k]) * 10 for k in ("x", "y", "z")] for p in positions]
        lines = topology.splitlines()
        atoms = [line for line in lines if line[:6] in {"ATOM  ", "HETATM"}]
        if not atoms or len(xyz) != len(atoms) or sum(line.startswith("MODEL ") for line in lines) > 1:
            return None, "The State positions do not match the retained topology's atom count."
        if not all(math.isfinite(v) for point in xyz for v in point):
            return None, "The retained State has non-finite coordinates."
        formatted = [[f"{v:8.3f}" for v in point] for point in xyz]
        if any(len(value) != 8 for point in formatted for value in point):
            return None, "The State coordinates cannot be represented safely as PDB."
        box = state.find("PeriodicBoxVectors")
        if box is None:
            return None, "The retained State has no readable periodic box vectors."
        vectors = [[float(box.find(name).attrib[k]) * 10 for k in ("x", "y", "z")]
                   for name in ("A", "B", "C")]
        if not all(math.isfinite(v) for row in vectors for v in row) or any(
                vectors[i][i] <= 0 for i in range(3)) or any(
                abs(vectors[i][j]) > 1e-7 for i, j in ((0, 1), (0, 2), (1, 2))):
            return None, "The State's box vectors cannot be represented safely as PDB."
        lengths = [math.sqrt(sum(v*v for v in row)) for row in vectors]
        def angle(i: int, j: int) -> float:
            dot = sum(a*b for a, b in zip(vectors[i], vectors[j])) / (lengths[i]*lengths[j])
            return math.degrees(math.acos(max(-1.0, min(1.0, dot))))
        cell = [f"{v:9.3f}" for v in lengths] + [f"{angle(i,j):7.2f}" for i,j in ((1,2),(0,2),(0,1))]
        if any(len(v) != (9 if i < 3 else 7) for i,v in enumerate(cell)):
            return None, "The State's cell cannot be represented safely as PDB."
        answer = ["CRYST1" + "".join(cell) + " P 1           1"]
        at = 0
        for line in lines:
            if line[:6] == "CRYST1":
                continue
            if line[:6] in {"ATOM  ", "HETATM"}:
                line = line[:30] + "".join(formatted[at]) + line[54:]
                at += 1
            answer.append(line)
        return "\n".join(answer) + "\n", ""
    except (ElementTree.ParseError, KeyError, ValueError, TypeError, AttributeError, OverflowError):
        return None, "The retained State or topology is malformed."


def stage_snapshot(root: str | Path, stage: str, version: str = "", *,
                   with_solvent: bool = False) -> tuple[bytes | None, dict[str, Any]]:
    """Read only a fixed stage candidate; caller text never names a path."""
    root = Path(root).resolve()
    try:
        with confined_records(root, lambda path: _public(root, path)) as scope:
            entry = next((e for e in _entries(root) if e["id"] == stage), None)
            if scope.denied:
                return None, {"reason": "A stage record could not be confined to the study's public files."}
            if entry is None or not entry["available"]:
                return None, {"reason": (entry or {}).get("reason") or "This stage has no retained coordinates."}
            if version and version != entry["version"]:
                return None, {"reason": "The stage files changed. Refresh the timeline before opening them."}
            pdb = _text(root / entry["source"])
            if pdb is None:
                return None, {"reason": "The retained stage source exceeds the display size limit."}
            if entry["topology"]:
                topology = _text(root / entry["topology"])
                if topology is None:
                    return None, {"reason": "The retained topology exceeds the display size limit."}
                pdb, reason = state_coordinates(topology, pdb)
                if pdb is None:
                    return None, {"reason": reason}
            if not any(line[:6] in {"ATOM  ", "HETATM"} for line in pdb.splitlines()):
                return None, {"reason": "The retained stage source contains no PDB atoms."}
            sources = [entry["source"]] + ([entry["topology"]] if entry["topology"] else [])
            if entry["version"] != _version(root, sources):
                return None, {"reason": "The stage files changed while being read. Try again."}
            if not with_solvent:
                from fastmdxplora.gui.live_frames import dashboard_display_pdb

                pdb = dashboard_display_pdb(pdb)
                if not any(line[:6] in {"ATOM  ", "HETATM"} for line in pdb.splitlines()):
                    return None, {"reason": "Only solvent coordinates remain. Include water and ions to inspect this stage."}
            return pdb.encode("utf-8"), entry
    except (OSError, UnicodeError, OutsideWorkspace):
        return None, {"reason": "Stage coordinates could not be read within the study's public files."}
