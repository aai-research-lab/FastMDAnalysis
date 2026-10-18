"""Inspect preparation observations without repeating a preparation operation."""
from __future__ import annotations

import hashlib
import json
import math
import re
from pathlib import Path
from typing import Any

from fastmdxplora.record_access import confined_records, open_text, stat
from fastmdxplora.refusals import OutsideWorkspace


AUDIT = "setup/preparation_audit.json"
MAX_AUDIT = 2_000_000
MAX_SNAPSHOT = 16 * 1024 * 1024


def _public(root: Path, path: Path) -> bool:
    from fastmdxplora.gui.server import _private

    try:
        return not _private(path.relative_to(root).parts)
    except ValueError:
        return False


def _read(path: Path, limit: int) -> str:
    if stat(path).st_size > limit:
        return ""
    with open_text(path, encoding="utf-8", newline="") as stream:
        text = stream.read(limit + 1)
    return text if len(text) <= limit else ""


def _audit(root: Path) -> tuple[dict[str, Any], str]:
    raw = _read(root / AUDIT, MAX_AUDIT)
    value = json.loads(raw)
    if not isinstance(value, dict) or value.get("version") != 1:
        return {}, ""
    if not isinstance(value.get("events"), list) or not isinstance(value.get("sources"), dict):
        return {}, ""
    if len(value["events"]) > 500 or len(value["sources"]) > 1000:
        return {}, ""
    return value, hashlib.sha256(raw.encode("utf-8")).hexdigest()


def _event(row: Any) -> dict[str, Any] | None:
    if not isinstance(row, dict) or not re.fullmatch(r"event-\d{4}", str(row.get("id", ""))):
        return None
    details = row.get("details")
    return {"id": row["id"], "operation": str(row.get("operation") or "Recorded operation")[:200],
            "before": row.get("before"), "after": row.get("after"),
            "details": details if isinstance(details, dict) else {},
            "reason": str(row.get("reason") or "")[:4000]}


def preparation_payload(root: str | Path) -> dict[str, Any]:
    """A bounded list of observations, kept distinct from requested choices."""
    root = Path(root).resolve()
    try:
        with confined_records(root, lambda path: _public(root, path)) as scope:
            audit, version = _audit(root)
            events = [event for row in audit.get("events", []) if (event := _event(row))]
        if scope.denied or not audit:
            return {"available": False, "reason": "No supported preparation audit could be read within this study."}
        warnings = audit.get("warnings")
        warnings = warnings if isinstance(warnings, list) else []
        return {"available": True, "version": version, "status": str(audit.get("status") or "unknown"),
                "notice": str(audit.get("notice") or "Only recorded preparation observations are shown."),
                "warnings": [str(s)[:4000] for s in warnings[:30]],
                "events": events}
    except (OSError, ValueError, UnicodeError, OutsideWorkspace):
        return {"available": False, "reason": "No supported preparation audit could be read within this study."}


def pdb_atoms(pdb: str) -> tuple[list[dict[str, Any]], str]:
    """Unique author identities for comparison, independent of atom serials.

    Models and alternate conformations need an explicit correspondence that
    the audit does not record, so they are never guessed from line numbers.
    """
    lines = pdb.splitlines()
    if sum(line.startswith("MODEL ") for line in lines) > 1:
        return [], "Several PDB models cannot be compared by one atom identity."
    atoms = []
    seen = set()
    try:
        for line in lines:
            if line[:6] not in {"ATOM  ", "HETATM"}:
                continue
            if len(line) < 54 or line[16:17].strip():
                return [], "Alternate or incomplete PDB atom records cannot be compared safely."
            number = line[22:26].strip()
            if not re.fullmatch(r"-?\d+", number):
                return [], "A PDB residue number cannot be mapped to its author identity."
            identity = (line[21:22].strip(), str(int(number)), line[26:27].strip(),
                        line[17:20].strip(), line[12:16].strip())
            if not identity[3] or not identity[4] or identity in seen:
                return [], "PDB atom identities are missing or ambiguous."
            point = [float(line[a:b]) for a, b in ((30, 38), (38, 46), (46, 54))]
            if not all(math.isfinite(v) for v in point):
                return [], "The retained PDB has non-finite coordinates."
            seen.add(identity)
            atoms.append({"identity": list(identity), "element": line[76:78].strip(), "xyz": point})
    except ValueError:
        return [], "The retained PDB has unreadable coordinates."
    return (atoms, "") if atoms else ([], "The retained snapshot contains no PDB atoms.")


def _snapshot(root: Path, sources: dict[str, Any], source_id: Any) -> tuple[dict[str, Any] | None, str]:
    if not isinstance(source_id, str) or not re.fullmatch(r"source-\d{4}", source_id):
        return None, "This observation has no retained structure."
    source = sources.get(source_id)
    if not isinstance(source, dict):
        return None, "The audit does not name this structure."
    if source.get("format") != "pdb" or source.get("snapshot") != f"setup/audit/{source_id}.pdb":
        return None, "A supported PDB snapshot was not retained for this observation."
    path = source["snapshot"]
    pdb = _read(root / path, MAX_SNAPSHOT)
    if not pdb or hashlib.sha256(pdb.encode("utf-8")).hexdigest() != source.get("sha256"):
        return None, "The retained snapshot is missing, changed or exceeds the display limit."
    atoms, reason = pdb_atoms(pdb)
    if reason:
        return None, reason
    return {"source": path, "pdb": pdb, "atoms": atoms}, ""


def preparation_change(root: str | Path, event_id: str, version: str) -> dict[str, Any]:
    """Check saved source hashes and compare identities for one observation."""
    root = Path(root).resolve()
    try:
        with confined_records(root, lambda path: _public(root, path)) as scope:
            audit, current = _audit(root)
            if not current or current != version:
                return {"available": False, "reason": "The preparation record changed. Read it again before inspecting a change."}
            events = [event for row in audit["events"] if (event := _event(row))]
            found = [event for event in events if event["id"] == event_id]
            if len(found) != 1:
                return {"available": False, "reason": "The preparation record does not identify one such observation."}
            event = found[0]
            before, before_reason = _snapshot(root, audit["sources"], event["before"])
            after, after_reason = _snapshot(root, audit["sources"], event["after"])
            if scope.denied:
                return {"available": False, "reason": "A preparation snapshot is outside this study's public files."}
        answer = {"available": bool(before or after), "event": event, "version": current,
                  "before": before, "after": after, "before_reason": before_reason,
                  "after_reason": after_reason, "compared": bool(before and after)}
        if before and after:
            old = {tuple(atom["identity"]) for atom in before["atoms"]}
            new = {tuple(atom["identity"]) for atom in after["atoms"]}
            answer["added"] = [atom["identity"] for atom in after["atoms"] if tuple(atom["identity"]) not in old]
            answer["removed"] = [atom["identity"] for atom in before["atoms"] if tuple(atom["identity"]) not in new]
            answer["unchanged"] = len(old & new)
        return answer
    except (OSError, ValueError, UnicodeError, OutsideWorkspace):
        return {"available": False, "reason": "The recorded preparation change could not be read within this study."}
