"""Problems already recorded by a study, without running or changing it."""

from __future__ import annotations

import json
import logging
import math
import re
from pathlib import Path
from typing import Any, Callable

from fastmdxplora.gui.analysis_overview import overview_of
from fastmdxplora.gui.series import NAME
from fastmdxplora.gui.telemetry import (
    EVENTS_FILE,
    METRICS_FILE,
    analyze_health,
    read_events,
    read_metrics,
    status_as_it_stands,
)
from fastmdxplora.gui.workspace import study_key
from fastmdxplora.record_access import confined_records, is_dir, is_file, read_text
from fastmdxplora.refusals import code_for
from fastmdxplora.remedies import remedy_for


_LOG = logging.getLogger(__name__)

#: Said under "What would fix it" where the records hold no remedy.
NO_REMEDY = "No remedy is recorded."

#: The run's narration, where its warnings are written (`TelemetryWriter.event`).
EVENTS = f"simulation/{EVENTS_FILE}"

#: Every line of the narration: each recorded warning is listed, not the
#: newest hundred lines' worth.
_EVERY_EVENT = 1_000_000

#: The warning the orchestrator writes beside a stop (`orchestrator.py`),
#: already listed as the stop itself.
_PHASE_STOPPED = re.compile(r"\w+ phase stopped")


def admitted_for_problems(root: Path) -> Callable[[Path], bool]:
    """The study's public files, and the run's own process record.

    The process record (`RUN_PROCESS_FILE`) is hidden, so the public rule
    refuses it, and without it a run whose record says it is going cannot
    be told from one whose process is gone. Only that file is admitted, in
    the study or in one of its pieces (`segment-*`), where
    `how_the_run_ended` reads it.
    """
    from fastmdxplora.gui.server import _private, _public_records
    from fastmdxplora.orchestrator import RUN_PROCESS_FILE

    public = _public_records(root)

    def admitted(path: Path) -> bool:
        try:
            parts = path.relative_to(root).parts
        except ValueError:
            return False
        if parts and parts[-1] == RUN_PROCESS_FILE and (
                len(parts) == 1
                or (len(parts) == 2 and parts[0].startswith("segment-")
                    and not _private(parts[:1]))):
            return True
        return public(path)

    return admitted


def problems_payload(root: str | Path, *,
                     artifacts: list[dict[str, Any]] | None = None,
                     admit_evidence: bool = False) -> dict[str, Any]:
    """Read existing records within the open study's public file boundary.

    Explicit ``artifacts`` further restricts evidence controls. The route
    may instead admit just the candidates used by this panel, without
    recursively listing the study. Missing records remain unknown, and
    ``read`` and ``unavailable`` name the records read and those that
    could not be.
    """
    root = Path(root).resolve()
    key = study_key(root)
    if not root.is_dir():
        return {"ok": True, "study": key, "state": "gone", "items": [], "available": False,
                "read": [], "unavailable": [],
                "notes": ["The study's folder is no longer there, so its records "
                          "cannot be read."]}
    try:
        with confined_records(root, admitted_for_problems(root)) as scope:
            answer = _read_problems(root, artifacts=artifacts, admit_evidence=admit_evidence)
    except OSError:
        _LOG.debug("Could not open the study's folder for reading", exc_info=True)
        return {"ok": True, "study": key, "state": "unreadable", "items": [],
                "available": False, "read": [], "unavailable": [],
                "notes": ["The study's folder could not be opened for reading, so its "
                          "problems are unknown."]}
    if scope.denied:
        # Even if a legacy reader caught the I/O refusal and used a
        # fallback, never present that fallback as a checked aggregate.
        return {"ok": True, "study": key, "items": [], "available": False,
                "read": [], "unavailable": sorted(scope.denied),
                "notes": ["Problems were not read because a source or link could not be "
                          "confined to the study's allowed files: "
                          + ", ".join(sorted(scope.denied)) + "."]}
    return answer


def _number(value: Any) -> float | None:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return None
    return float(value) if math.isfinite(value) else None


def _more_sampling(shortfall: Any) -> str:
    """How much more a withheld mean asked for, in its own record's numbers
    (`statistics.Shortfall.as_record`); empty where it recorded none."""
    if not isinstance(shortfall, dict):
        return ""
    from fastmdxplora.simulation.sampling_ask import length_said

    target = _number(shortfall.get("target_independent"))
    frames = _number(shortfall.get("more_frames"))
    more_ns = _number(shortfall.get("more_ns"))
    if target is None or target <= 0 or frames is None or frames <= 0:
        return ""
    floor = bool(shortfall.get("lower_bound"))
    counted = f"{int(frames)} further frame{'s' if int(frames) != 1 else ''}"
    amount = (f"{length_said(more_ns)} more production ({counted})"
              if more_ns is not None and more_ns > 0 else counted)
    said = (f"{'At least ' if floor else ''}{amount} should give "
            f"{target:g} independent samples.")
    if floor:
        said += (" This is a floor: the run is too short to resolve its own "
                 "correlation time, so run that much and estimate again.")
    return said[:1].upper() + said[1:]


def _quantity_record(findings: dict[str, Any], key: str) -> dict[str, Any] | None:
    """The record a quantity of the Overview was read from (`overview_of`)."""
    found = findings.get(key)
    if found is None and "." in key:
        group, inner = key.split(".", 1)
        found = findings.get(group, {}).get(inner) if isinstance(
            findings.get(group), dict) else None
    return found if isinstance(found, dict) else None


def _read_problems(root: Path, *, artifacts: list[dict[str, Any]] | None,
                   admit_evidence: bool) -> dict[str, Any]:
    notes: list[str] = []
    items: list[dict[str, str]] = []
    read: list[str] = []
    unavailable: list[str] = []
    allowed = {entry.get("path") for entry in artifacts or []
               if isinstance(entry, dict) and isinstance(entry.get("path"), str)}
    available = False

    def add(identifier: str, title: str, reason: str, fix: str, *,
            code: str = "", files: tuple[str, ...] = ()) -> None:
        item = {"id": identifier, "title": title, "reason": reason, "fix": fix or NO_REMEDY}
        if code:
            item["code"] = code
        for path in files:
            if (path in allowed or (artifacts is None and admit_evidence)) and is_file(root / path):
                item["file"] = path
                break
        items.append(item)

    def unknown(relative: str, note: str) -> None:
        if relative not in unavailable:
            unavailable.append(relative)
        notes.append(note)

    def record(relative: str, missing: str = "") -> dict[str, Any] | None:
        path = root / relative
        try:
            value = json.loads(read_text(path, encoding="utf-8"))
        except FileNotFoundError:
            if missing:
                unknown(relative, f"{relative} is not there; {missing}")
            return None
        except (OSError, ValueError):
            unknown(relative, f"{relative} could not be read; its problems are unknown.")
            return None
        if not isinstance(value, dict) or not value:
            unknown(relative, f"{relative} has no usable record; its problems are unknown.")
            return None
        if relative not in read:
            read.append(relative)
        return value

    manifest = record("manifest.json", "it is written when a run ends, so refusals and "
                                       "phase results are unknown.")
    phases = manifest.get("phases") if manifest else None
    if manifest is not None and (not isinstance(phases, list) or not phases):
        unknown("manifest.json", "The manifest has no readable phase records.")
    for index, phase in enumerate(phases if isinstance(phases, list) else []):
        if not isinstance(phase, dict) or not phase.get("status"):
            unknown("manifest.json", "A manifest phase has no readable status.")
            continue
        available = True
        refusal = phase.get("refusal")
        failed = str(phase["status"]).lower() in {"error", "failed", "stopped", "interrupted"}
        if not failed and not refusal:
            continue
        name = str(phase.get("name") or phase.get("phase") or "Study")
        reason = str(phase.get("message") or "No reason was recorded.")
        fix = NO_REMEDY
        code = ""
        if isinstance(refusal, dict):
            reason = str(refusal.get("message") or phase.get("message") or reason)
            code = str(refusal.get("code") or "")
            try:
                # A retry depends on whether the run wrote a checkpoint, so
                # its remedy is asked of the study (read, never priced).
                # Otherwise the registry alone answers, without a sampling
                # estimate or a Config built from the study.
                retry = code_for(code).retryable
                fix = (remedy_for(refusal, study=root, price=False) if retry
                       else remedy_for(refusal)).fix
            except Exception:  # noqa: BLE001 - retain the refusal when its remedy cannot be read
                _LOG.debug("Could not read a recorded phase remedy", exc_info=True)
                notes.append(f"The remedy for {name} could not be read.")
        elif refusal is not None:
            notes.append(f"The refusal record for {name} could not be read.")
        state = (code.rsplit(".", 1)[-1] if code in
                 {"simulation.run.stopped", "simulation.run.interrupted"} else phase["status"])
        add(f"phase:{index}", f"{name.capitalize()}: {state}", reason, fix,
            code=code, files=("manifest.json",))

    raw_status = record("simulation/live_status.json", "the run's health is unknown.")
    try:
        # As the Overview and the sidebar say it: a run whose record says it
        # is going but whose process is gone, or silent for a day, ended.
        status = status_as_it_stands(root, include_campaign=False)
        metrics = read_metrics(root)
        if metrics:
            read.append(f"simulation/{METRICS_FILE}")
        health = analyze_health(status, metrics)
        state = health.get("state")
        if status.get("status"):
            available = True
        elif raw_status is not None:
            unknown("simulation/live_status.json", "The telemetry record has no run status.")
        if state in {"warning", "failed", "stopped", "interrupted"}:
            reason = str(health.get("message") or "No reason was recorded.")
            fix = str(health.get("explanation") or NO_REMEDY)
            if state in {"stopped", "interrupted"}:
                if status.get("stage") in {"analysis", "report"}:
                    fix = "Use Analyze again or Write it again on the Overview for the phase that stopped."
                else:
                    # Checkpoint-aware wording, without pricing an action
                    # that the read-only panel never offers to run.
                    fix = remedy_for({"code": f"simulation.run.{state}", "message": reason},
                                     study=root, price=False).fix
            duplicate = next((item for item in items if
                item["reason"].strip() == reason.strip()
                or (state in {"stopped", "interrupted"}
                    and item.get("code") == f"simulation.run.{state}")), None)
            if duplicate is not None and state in {"stopped", "interrupted"}:
                duplicate["fix"] = fix
            if duplicate is None:
                add(f"health:{state}", str(health.get("headline") or f"Run {state}"),
                    reason, fix,
                    files=("simulation/live_status.json", f"simulation/{METRICS_FILE}"))
        elif state == "unknown" and raw_status is not None:
            notes.append(str(health.get("message") or "Run health is unavailable."))
    except Exception:  # noqa: BLE001 - partial records must not hide the other problems
        _LOG.debug("Could not read recorded run health", exc_info=True)
        unknown("simulation/live_status.json",
                "Run health could not be read; its problems are unknown.")

    try:
        # Each warning the run wrote in its narration (runner.py and
        # orchestrator.py write them with `TelemetryWriter.event`), once.
        events = read_events(root, limit=_EVERY_EVENT)
        if is_file(root / EVENTS):
            read.append(EVENTS)
        stop = any(item["id"] == "health:stopped" or item.get("code") == "simulation.run.stopped"
                   for item in items)
        seen: set[str] = set()
        for event in events:
            message = str(event.get("message") or "").strip()
            if str(event.get("level") or "").lower() != "warning" or not message:
                continue
            if message in seen:
                continue
            seen.add(message)
            if any(item["reason"].strip() == message for item in items):
                continue
            if stop and _PHASE_STOPPED.fullmatch(message):
                continue
            add(f"warning:{len(seen)}", "Recorded run warning", message, NO_REMEDY,
                files=(EVENTS,))
    except Exception:  # noqa: BLE001 - an unreadable narration must not hide the rest
        _LOG.debug("Could not read the run's recorded warnings", exc_info=True)
        unknown(EVENTS, "The run's recorded warnings could not be read; they are unknown.")

    analysis_expected = is_dir(root / "analysis") or any(
        isinstance(phase, dict) and str(phase.get("name") or "") == "analysis"
        for phase in (phases if isinstance(phases, list) else []))
    analysis = record("analysis/analysis_manifest.json",
                      "a failed analysis would not be known." if analysis_expected else "")
    results = analysis.get("results") if analysis else None
    invalid_names = isinstance(results, dict) and any(
        not NAME.fullmatch(name) for name in results)
    if analysis is not None and (not isinstance(results, dict) or invalid_names):
        unknown("analysis/analysis_manifest.json",
                "The analysis manifest has unreadable result entries.")
    elif isinstance(results, dict) and any(not isinstance(result, dict)
                                          or not result.get("status")
                                          for result in results.values()):
        unknown("analysis/analysis_manifest.json",
                "The analysis manifest has unreadable result entries.")
    try:
        # Do not follow an analysis name outside its folder, even in an old
        # or hand-edited manifest. The existing reader assumes valid names.
        overview = {"rows": []} if invalid_names else overview_of(root, include_time=False, include_data=False)
        rows = overview.get("rows") or []
        if not rows:
            notes.append("No analysis results are available to check.")
        for row in rows:
            name, title = row["analysis"], row["title"]
            options = f"analysis/{name}/options.json"
            found = record(options)
            if found is not None and not isinstance(found.get("findings"), dict):
                unknown(options, f"{title} has no readable findings record.")
                found = None
            if found is not None and "mean" in found["findings"] and not isinstance(
                    found["findings"]["mean"], dict):
                unknown(options, f"{title} has an unreadable mean record in {options}; "
                                 "its problems are unknown.")
            if found is None and row.get("status") not in {"failed", "skipped", "running"}:
                unknown(options, f"{title} has no readable findings record.")
            available = available or found is not None or row.get("status") == "failed"
            if row.get("status") == "failed":
                add(f"analysis:{name}:failed", f"{title}: failed",
                    str(row.get("message") or "No reason was recorded."),
                    NO_REMEDY, files=("analysis/analysis_manifest.json", options))
            for quantity in row.get("quantities") or []:
                earlier = bool(quantity.get("judged_earlier"))
                if quantity.get("determined") and not earlier:
                    continue
                label = quantity.get("label") or title
                named = title if label == title else f"{title} / {label}"
                if earlier:
                    add(f"analysis:{name}:{quantity['key']}", f"{named}: earlier rules",
                        "This mean was judged by earlier rules; its recorded verdict "
                        "has not been checked under the current rules.",
                        "Analyze again to apply the current rules.", files=(options,))
                    continue
                why = str(quantity.get("why") or "No reason was recorded.")
                own = _quantity_record(found["findings"], str(quantity["key"])) \
                    if found is not None else None
                fix = ""
                if own is not None and str(own.get("not_a_measurement") or "") == why:
                    # The record's own reason, so its own remedy: how much
                    # more sampling it asked for, where it asked.
                    fix = _more_sampling(own.get("shortfall"))
                # The same in every frame is said as the Overview says it,
                # not as a value missing.
                verdict = (str(quantity.get("said")) if quantity.get("constant")
                           and quantity.get("said") else "not determined")
                add(f"analysis:{name}:{quantity['key']}", f"{named}: {verdict}",
                    why, fix, files=(options,))
    except Exception:  # noqa: BLE001 - an unreadable analysis must not erase run problems
        _LOG.debug("Could not read recorded analysis results", exc_info=True)
        unknown("analysis", "Some analysis results could not be read; their problems "
                            "are unknown.")

    if is_file(root / "batch_manifest.json"):
        notes.append("This list covers the open folder's records. Open each campaign run "
                     "to check its health, refusals and analysis results.")
    if not available:
        notes.append("No usable run or analysis records are available.")
    return {"ok": True, "study": study_key(root), "items": items,
            "notes": list(dict.fromkeys(notes)), "available": available,
            "read": list(dict.fromkeys(read)), "unavailable": unavailable}
