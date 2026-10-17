"""Problems already recorded by a study, without running or changing it."""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Any

from fastmdxplora.gui.analysis_overview import overview_of
from fastmdxplora.gui.series import NAME
from fastmdxplora.gui.telemetry import analyze_health, read_metrics, status_as_it_stands
from fastmdxplora.gui.workspace import study_key
from fastmdxplora.record_access import confined_records, is_file, read_text
from fastmdxplora.refusals import code_for
from fastmdxplora.remedies import remedy_for


_LOG = logging.getLogger(__name__)


def problems_payload(root: str | Path, *,
                     artifacts: list[dict[str, Any]] | None = None,
                     admit_evidence: bool = False) -> dict[str, Any]:
    """Read existing records within the open study's public file boundary.

    Explicit ``artifacts`` further restricts evidence controls. The route
    may instead admit just the candidates used by this panel, without
    recursively listing the study. Missing records remain unknown.
    """
    from fastmdxplora.gui.server import _private

    root = Path(root).resolve()

    def public(path: Path) -> bool:
        try:
            return not _private(path.relative_to(root).parts)
        except ValueError:
            return False

    with confined_records(root, public) as scope:
        answer = _read_problems(root, artifacts=artifacts, admit_evidence=admit_evidence)
    if scope.denied:
        # Even if a legacy reader caught the I/O refusal and used a
        # fallback, never present that fallback as a checked aggregate.
        return {"ok": True, "study": study_key(root), "items": [], "available": False,
                "notes": ["Problems were not read because a source or link could not be "
                          "confined to the study's allowed files. Review the recorded "
                          "sources before checking again."]}
    return answer


def _read_problems(root: Path, *, artifacts: list[dict[str, Any]] | None,
                   admit_evidence: bool) -> dict[str, Any]:
    notes: list[str] = []
    items: list[dict[str, str]] = []
    allowed = {entry.get("path") for entry in artifacts or []
               if isinstance(entry, dict) and isinstance(entry.get("path"), str)}
    available = False

    def add(identifier: str, title: str, reason: str, fix: str, *,
            code: str = "", files: tuple[str, ...] = ()) -> None:
        item = {"id": identifier, "title": title, "reason": reason, "fix": fix}
        if code:
            item["code"] = code
        for path in files:
            if (path in allowed or (artifacts is None and admit_evidence)) and is_file(root / path):
                item["file"] = path
                break
        items.append(item)

    def record(relative: str) -> dict[str, Any] | None:
        path = root / relative
        try:
            value = json.loads(read_text(path, encoding="utf-8"))
        except FileNotFoundError:
            return None
        except (OSError, ValueError):
            notes.append(f"{relative} could not be read; its problems are unknown.")
            return None
        if not isinstance(value, dict) or not value:
            notes.append(f"{relative} has no usable record; its problems are unknown.")
            return None
        return value

    manifest = record("manifest.json")
    phases = manifest.get("phases") if manifest else None
    if manifest is not None and (not isinstance(phases, list) or not phases):
        notes.append("The manifest has no readable phase records.")
    for index, phase in enumerate(phases if isinstance(phases, list) else []):
        if not isinstance(phase, dict) or not phase.get("status"):
            notes.append("A manifest phase has no readable status.")
            continue
        available = True
        refusal = phase.get("refusal")
        failed = str(phase["status"]).lower() in {"error", "failed", "stopped", "interrupted"}
        if not failed and not refusal:
            continue
        name = str(phase.get("name") or phase.get("phase") or "Study")
        reason = str(phase.get("message") or "No reason was recorded.")
        fix = "Review the recorded reason and this phase's results."
        code = ""
        if isinstance(refusal, dict):
            reason = str(refusal.get("message") or phase.get("message") or reason)
            code = str(refusal.get("code") or "")
            try:
                # Without a study argument this reads the registry, not a new
                # sampling estimate or a priced command built from Config.
                remedy = remedy_for(refusal)
                # A retry depends on the checkpoint and whether another
                # process still runs. A registry alone cannot settle that.
                if not code_for(remedy.code).retryable:
                    fix = remedy.fix
            except Exception:  # noqa: BLE001 - retain the refusal when its remedy cannot be read
                _LOG.debug("Could not read a recorded phase remedy", exc_info=True)
                notes.append(f"The remedy for {name} could not be read.")
        elif refusal is not None:
            notes.append(f"The refusal record for {name} could not be read.")
        state = (code.rsplit(".", 1)[-1] if code in
                 {"simulation.run.stopped", "simulation.run.interrupted"} else phase["status"])
        add(f"phase:{index}", f"{name.capitalize()}: {state}", reason, fix,
            code=code, files=("manifest.json",))

    raw_status = record("simulation/live_status.json")
    try:
        status = status_as_it_stands(root, include_campaign=False)
        # Folder discovery belongs to the workspace. This panel states
        # missing telemetry directly rather than searching child folders.
        health = analyze_health(status, read_metrics(root))
        state = health.get("state")
        if status.get("status"):
            available = True
        elif raw_status is not None:
            notes.append("The telemetry record has no run status.")
        if state in {"warning", "failed", "stopped", "interrupted"}:
            reason = str(health.get("message") or "No reason was recorded.")
            fix = str(health.get("explanation") or
                      "Review the recorded reason and the run's results.")
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
                    files=("simulation/live_status.json", "simulation/live_metrics.csv"))
        elif state == "unknown":
            notes.append(str(health.get("message") or "Run health is unavailable."))
        for source in (raw_status or {}, status):
            warning = source.get("latest_warning")
            if not isinstance(warning, str) or not warning.strip():
                continue
            if not any(item["reason"].strip() == warning.strip() for item in items):
                add("health:recorded-warning", "Recorded warning", warning,
                    "Review the recorded warning and the run's results.",
                    files=("simulation/live_status.json",))
    except Exception:  # noqa: BLE001 - partial records must not hide the other problems
        _LOG.debug("Could not read recorded run health", exc_info=True)
        notes.append("Run health could not be read; its problems are unknown.")

    analysis = record("analysis/analysis_manifest.json")
    results = analysis.get("results") if analysis else None
    invalid_names = isinstance(results, dict) and any(
        not NAME.fullmatch(name) for name in results)
    if analysis is not None and (not isinstance(results, dict) or invalid_names):
        notes.append("The analysis manifest has unreadable result entries.")
    elif isinstance(results, dict) and any(not isinstance(result, dict)
                                          or not result.get("status")
                                          for result in results.values()):
        notes.append("The analysis manifest has unreadable result entries.")
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
                notes.append(f"{title} has no readable findings record.")
                found = None
            if found is not None and "mean" in found["findings"] and not isinstance(
                    found["findings"]["mean"], dict):
                notes.append(f"{title} has an unreadable mean record in {options}; "
                             "its problems are unknown.")
            if found is None and row.get("status") not in {"failed", "skipped", "running"}:
                notes.append(f"{title} has no readable findings record.")
            available = available or found is not None or row.get("status") == "failed"
            if row.get("status") == "failed":
                add(f"analysis:{name}:failed", f"{title}: failed",
                    str(row.get("message") or "No reason was recorded."),
                    "Review the recorded reason and this analysis's results.",
                    files=("analysis/analysis_manifest.json", options))
            for quantity in row.get("quantities") or []:
                earlier = bool(quantity.get("judged_earlier"))
                if quantity.get("determined") and not earlier:
                    continue
                label = quantity.get("label") or title
                reason = ("This mean was judged by earlier rules; its recorded verdict "
                          "has not been checked under the current rules." if earlier else
                          str(quantity.get("why") or "No reason was recorded."))
                named = title if label == title else f"{title} / {label}"
                add(f"analysis:{name}:{quantity['key']}",
                    f"{named}: {'earlier rules' if earlier else 'not determined'}",
                    reason, "Analyze again to apply the current rules." if earlier else
                    "Review the recorded reason and this analysis's results.", files=(options,))
    except Exception:  # noqa: BLE001 - an unreadable analysis must not erase run problems
        _LOG.debug("Could not read recorded analysis results", exc_info=True)
        notes.append("Some analysis results could not be read; their problems are unknown.")

    if is_file(root / "batch_manifest.json"):
        notes.append("This list covers the open folder's records. Open each campaign run "
                     "to check its health, refusals and analysis results.")
    if not available:
        notes.append("No usable run or analysis records are available.")
    return {"ok": True, "study": study_key(root), "items": items,
            "notes": list(dict.fromkeys(notes)), "available": available}
