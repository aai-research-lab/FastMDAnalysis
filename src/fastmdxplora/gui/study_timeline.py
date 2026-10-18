"""What happened to a study, newest first, from its own records.

A study's folder holds its history in pieces: the Manifest's record of each
phase as it last ran (when, how long, by which version on which computer),
the record and the output of a phase a run again replaced
(`previous/<phase>`, :mod:`fastmdxplora.replaced`), and the pieces an
extended study was carried on in (`segment-*`). This reads them into one
list, an editor's Timeline of a file: each step with when it happened, how
long it took, what made it, and what it changed, where the records say: the
means an analysis run again moved, the settings a piece was carried on
with. Nothing is said that the records do not hold: which person or Agent
asked for a step is not recorded, so it is not said.

Every record is read through `record_access`, so the route's read held to
the study (`server._confined`) refuses one linked out of it, and names it
as not read: a mean whose record was not read is never said added or
dropped, nor a piece's settings as the study's.
"""

from __future__ import annotations

import json
import math
import os
from datetime import datetime
from pathlib import Path
from typing import Any

from fastmdxplora import record_access

__all__ = ["timeline_of"]

#: What each phase is said as, done.
DONE = {"setup": "Set up", "simulation": "Simulated", "analysis": "Analysed",
        "report": "Report written"}

#: A mean is said to have moved past this fraction of itself (or of its
#: error, where it has one): a value read again to the last digit is the
#: same value.
MOVED = 1e-9


def timeline_of(root: str | Path) -> dict[str, Any]:
    """The study's steps, newest first."""
    base = Path(root)
    if not record_access.is_dir(base):
        return {"ok": False, "reason": "No study is open."}
    events: list[dict[str, Any]] = []
    manifest = _json(base / "manifest.json")
    phases = manifest.get("phases") if isinstance(manifest, dict) else None
    for record in phases if isinstance(phases, list) else []:
        if isinstance(record, dict) and record.get("name") in DONE:
            events.append(_phase(base, record, replaced=False))
    previous = base / "previous"
    ran_again = {e["phase"] for e in events}
    for phase in DONE:
        record = _json(previous / phase / "phase_record.json")
        if isinstance(record, dict) and record.get("name") == phase:
            event = _phase(base, record, replaced=True, again=phase in ran_again)
            events.append(event)
            if phase == "analysis":
                again = next((e for e in events if e["phase"] == "analysis" and not e["replaced"]),
                             None)
                if again is not None:
                    again["changes"] = _means_moved(previous / "analysis", base / "analysis")
                    again["again"] = True
    events += _pieces(base)
    events.sort(key=lambda e: e.get("at") or "", reverse=True)
    return {"ok": True, "events": events}


def _phase(base: Path, record: dict[str, Any], *, replaced: bool,
           again: bool = False) -> dict[str, Any]:
    name = str(record.get("name"))
    made = record.get("produced_by") if isinstance(record.get("produced_by"), dict) else {}
    status = str(record.get("status") or "")
    started, finished = record.get("started_at"), record.get("finished_at")
    said = DONE[name] if status in ("ok", "completed", "") else f"{name.capitalize()} {status}"
    event = {
        "kind": "phase", "phase": name, "said": said, "status": status,
        "at": finished or started, "started": started,
        "seconds": _seconds(started, finished),
        "version": str(made.get("version") or ""), "host": str(made.get("host") or ""),
        "replaced": replaced, "again": False, "changes": None,
    }
    if replaced:
        # Kept by a run again, which moved it to previous/: its output is
        # read there, never as the study's. Replaced where it ran again;
        # set aside where it was written from what ran again and has not
        # been written since.
        event["said"] = said + (", then replaced" if again
                                else ", then set aside: written from what was run again")
        event["where"] = f"previous/{name}"
    detail = _detail(base / "previous" / name if replaced else base / name, name)
    if detail:
        event["detail"] = detail
    return event


def _detail(folder: Path, phase: str) -> str:
    """What a phase's record says it covered, in a few words."""
    if phase == "analysis":
        plan = _json(folder / "analysis_manifest.json")
        plan = plan.get("plan") if isinstance(plan, dict) else None
        if isinstance(plan, list) and plan:
            return f"{len(plan)} analys{'is' if len(plan) == 1 else 'es'}: " + ", ".join(
                map(str, plan[:6])) + (", ..." if len(plan) > 6 else "")
    return ""


def _means_moved(before: Path, after: Path) -> dict[str, Any]:
    """The means an analysis run again changed: each analysis in both whose
    mean moved, or became determined or not; those added and those dropped.
    One whose record was not read on either side is said not read, never
    added, dropped or kept."""
    from fastmdxplora.gui.report_dashboard import unit_of

    (old, unread_before), (new, unread_after) = _means(before), _means(after)
    not_read = sorted(unread_before | unread_after)
    for name in not_read:
        old.pop(name, None)
        new.pop(name, None)
    moved = []
    for name in sorted(set(old) & set(new)):
        a, b = old[name], new[name]
        value_a, value_b = _finite(a.get("mean")), _finite(b.get("mean"))
        error_a, error_b = _finite(a.get("standard_error")), _finite(b.get("standard_error"))
        determined_a = error_a is not None and not a.get("not_a_measurement")
        determined_b = error_b is not None and not b.get("not_a_measurement")
        scale = max(abs(value_a or 0.0), abs(value_b or 0.0), error_a or 0.0, error_b or 0.0, 1e-300)
        changed = (value_a is None) != (value_b is None) or (
            value_a is not None and abs(value_b - value_a) > MOVED * scale)
        if changed or determined_a != determined_b:
            moved.append({"analysis": name, "unit": unit_of(name, b),
                          "before": {"mean": value_a, "error": error_a if determined_a else None},
                          "after": {"mean": value_b, "error": error_b if determined_b else None}})
    return {"moved": moved, "added": sorted(set(new) - set(old)),
            "dropped": sorted(set(old) - set(new)),
            "kept": len(set(old) & set(new)) - len(moved), "not_read": not_read}


def _means(folder: Path) -> tuple[dict[str, dict[str, Any]], set[str]]:
    """Each analysis's mean over the trajectory, as its record keeps it; and
    the analyses whose record is there but was not read."""
    found: dict[str, dict[str, Any]] = {}
    not_read: set[str] = set()
    try:
        folders = [f for f in record_access.iterdir(folder) if record_access.is_dir(f)]
    except OSError:
        return found, not_read
    for each in folders:
        path = each / "options.json"
        if os.path.lexists(path) and not record_access.can_read(path):
            not_read.add(each.name)
    try:
        records = sorted(record_access.glob(folder, "*/options.json"))
    except OSError:
        return found, not_read
    for path in records:
        record = _json(path)
        findings = record.get("findings") if isinstance(record, dict) else None
        mean = findings.get("mean") if isinstance(findings, dict) else None
        if isinstance(mean, dict):
            found[path.parent.name] = mean
    return found, not_read


def _pieces(base: Path) -> list[dict[str, Any]]:
    """Each piece an extended study was carried on in, with the settings it
    was carried on with that differ from the study's own."""
    import re

    events = []
    study = _config(base)
    for folder in sorted(record_access.glob(base, "segment-*")):
        number = re.fullmatch(r"segment-(\d+)", folder.name)
        if not number or not record_access.is_dir(folder):
            continue
        manifest = _json(folder / "manifest.json")
        phases = manifest.get("phases") if isinstance(manifest, dict) else None
        simulated = next((p for p in phases or [] if isinstance(p, dict)
                          and p.get("name") == "simulation"), None)
        if simulated is None:
            continue
        made = simulated.get("produced_by") if isinstance(simulated.get("produced_by"), dict) else {}
        started, finished = simulated.get("started_at"), simulated.get("finished_at")
        settings = []
        piece = _config(folder)
        # Its settings, or the study's, not read: said so, never as a piece
        # carried on with the study's own.
        settings_not_read = any(
            os.path.lexists(where / "resolved_config.yml")
            and not record_access.can_read(where / "resolved_config.yml")
            for where in (base, folder))
        if study is not None and piece is not None:
            from fastmdxplora.config.diff import differences

            try:
                settings = [_said_in_its_unit(d.as_record()) for d in differences(study, piece)
                            if not d.where_only and not d.unrecorded
                            and d.setting not in _CARRYING][:12]
            except Exception:  # noqa: BLE001 - the step is said without them
                settings = []
        events.append({
            "kind": "piece", "phase": "simulation",
            "said": f"Carried on, piece {int(number.group(1))}",
            "status": str(simulated.get("status") or ""),
            "at": finished or started, "started": started,
            "seconds": _seconds(started, finished),
            "version": str(made.get("version") or ""), "host": str(made.get("host") or ""),
            "replaced": False, "again": False, "changes": None, "settings": settings,
            "settings_not_read": settings_not_read, "where": folder.name,
        })
    return events


#: What carrying a study on writes into its piece's settings, the same for
#: every piece: where it was prepared and resumes from, no minimisation or
#: equilibration again, and only the simulation run. Not what the piece was
#: carried on with.
_CARRYING = frozenset({
    "include_phase", "simulation.setup_from", "simulation.prepared_from",
    "simulation.resume_from", "simulation.resume_unsealed", "simulation.minimize",
    "simulation.nvt_steps", "simulation.npt_steps", "simulation.nvt_duration_ns",
    "simulation.npt_duration_ns",
})


def _said_in_its_unit(record: dict[str, Any]) -> dict[str, Any]:
    """A setting in nanoseconds also said as the pages say a length, in
    picoseconds under one nanosecond (``first_said``, ``second_said``)."""
    from fastmdxplora.simulation.sampling_ask import length_said

    if str(record.get("setting", "")).endswith("_ns"):
        for side in ("first", "second"):
            value = record.get(side)
            if isinstance(value, (int, float)) and not isinstance(value, bool):
                record[f"{side}_said"] = length_said(value)
    return record


def _config(folder: Path) -> dict[str, Any] | None:
    from fastmdxplora.config.diff import config_of

    try:
        config, _ = config_of(folder)
    except Exception:  # noqa: BLE001 - none recorded
        return None
    return config if isinstance(config, dict) else None


def _seconds(started: Any, finished: Any) -> float | None:
    try:
        a = datetime.fromisoformat(str(started))
        b = datetime.fromisoformat(str(finished))
    except (TypeError, ValueError):
        return None
    seconds = (b - a).total_seconds()
    return seconds if seconds >= 0 else None


def _finite(value: Any) -> float | None:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return None
    return float(value) if math.isfinite(value) else None


def _json(path: Path) -> Any:
    try:
        return json.loads(record_access.read_text(path, encoding="utf-8"))
    except (OSError, ValueError):
        return None
