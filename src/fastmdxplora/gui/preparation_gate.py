"""A local GUI handoff from audited setup to an explicitly accepted simulation."""
from __future__ import annotations

from copy import deepcopy
import hashlib
import json
from pathlib import Path
import shutil
import sys
from typing import Any
import uuid

from fastmdxplora.record_access import confined_records, open_text, stat
from fastmdxplora.refusals import OutsideWorkspace, StudyError
from fastmdxplora.gui.preparation_review import AUDIT, _audit, _public, _read, _snapshot


PLAN = "simulation_plan.yml"
RECORD = ".preparation-review.json"
CORE = ("manifest.json", "resolved_config.yml", "setup/setup_parameters.json",
        "setup/system.xml", "setup/state.xml", "setup/topology.pdb", AUDIT, PLAN)
MAX_CORE = 256 * 1024 * 1024


def scope_reason(config: dict[str, Any], *, hosted: bool = False) -> str:
    """Refuse unsupported workflows before writing or starting anything."""
    from fastmdxplora.gui.config_builder import phases_it_runs

    if hosted:
        return "Review after setup is available for local runs on this machine."
    systems = config.get("systems")
    if not isinstance(systems, list) or len(systems) != 1 or config.get("sweep"):
        return "Review after setup needs one system without a sweep."
    if not {"setup", "simulation"}.issubset(phases_it_runs(config)):
        return "Review after setup needs both Setup and Simulation in the Config."
    simulation = {**(config.get("simulation") or {}), **(systems[0].get("simulation") or {})}
    if any(simulation.get(key) for key in ("setup_from", "prepared_from", "resume_from", "extra_ns")):
        return "Review after setup prepares a new system; reused setups and continuations need their existing workflow."
    if config.get("budget_hours") is not None:
        return "Review after setup cannot split a study with a compute budget. Remove the budget or use the existing staged workflow."
    return ""


def _write(root: Path, record: dict[str, Any]) -> None:
    temporary = root / (RECORD + "." + uuid.uuid4().hex)
    try:
        with temporary.open("x", encoding="utf-8") as stream:
            json.dump(record, stream, indent=2)
        temporary.replace(root / RECORD)
    finally:
        temporary.unlink(missing_ok=True)


def arm(root: Path, config: dict[str, Any]) -> dict[str, Any]:
    """Keep the full study plan and return the Config for setup alone."""
    import yaml
    from fastmdxplora.gui.config_builder import export_config

    reason = scope_reason(config)
    if reason:
        raise StudyError(reason, code="config.option.inapplicable")
    exported = export_config(deepcopy(config), full=True)
    if not exported.get("ok"):
        raise StudyError(exported.get("error") or "The simulation plan is invalid.", code="config.option.inapplicable")
    plan = exported["yaml"]
    (root / PLAN).write_bytes(plan.encode("utf-8"))
    _write(root, {"version": 1, "status": "waiting", "plan_sha256": hashlib.sha256(plan.encode()).hexdigest()})
    setup = yaml.safe_load(plan)
    setup["include_phase"] = ["setup"]
    setup.pop("exclude_phase", None)
    setup.setdefault("setup", {})["preparation_audit"] = True
    # A system's own value takes precedence over the phase value.
    for system in setup["systems"]:
        system.setdefault("setup", {})["preparation_audit"] = True
    return setup


def _hash(path: Path) -> str:
    before = stat(path)
    if before.st_size > MAX_CORE:
        raise StudyError("A preparation source exceeds the review size limit.", code="setup.review.unavailable")
    digest = hashlib.sha256()
    size = 0
    with open_text(path, encoding="utf-8", newline="") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), ""):
            data = chunk.encode("utf-8")
            size += len(data)
            if size > MAX_CORE:
                raise StudyError("A preparation source exceeds the review size limit.", code="setup.review.unavailable")
            digest.update(data)
    after = stat(path)
    if (before.st_mtime_ns, before.st_ctime_ns, before.st_size) != (
            after.st_mtime_ns, after.st_ctime_ns, after.st_size):
        raise StudyError("The preparation changed while being read. Refresh the review.", code="setup.review.unavailable")
    return digest.hexdigest()


def review_state(root: str | Path, *, running: bool = False, sources: bool = False) -> dict[str, Any]:
    """A review version binds the displayed evidence to the retained plan."""
    root = Path(root).resolve()
    if not (root / RECORD).exists():
        return {"requested": False, "can_accept": False}
    if running:
        return {"requested": True, "can_accept": False, "status": "waiting",
                "reason": "Setup is still running. Simulation will wait for your acceptance."}
    try:
        with confined_records(root, lambda path: path == root / RECORD or _public(root, path)) as scope:
            record = json.loads(_read(root / RECORD, 64 * 1024))
            if not isinstance(record, dict) or record.get("version") != 1:
                raise StudyError("The preparation review record is unsupported.", code="setup.review.unavailable")
            if record.get("status") != "waiting":
                return {"requested": True, "can_accept": False, "status": record.get("status"),
                        "reason": ("Setup accepted. Simulation was started on this prepared study."
                                   if record.get("status") == "accepted" else
                                   "The simulation launch is being resolved. Inspect the current job before trying again.")}
            hashes = {name: _hash(root / name) for name in CORE}
            if hashes[PLAN] != record.get("plan_sha256"):
                raise StudyError("The retained simulation plan changed. Start a new reviewed study for a different plan.", code="setup.review.changed")
            manifest = json.loads(_read(root / "manifest.json", 2_000_000))
            phases = manifest.get("phases") if isinstance(manifest, dict) else None
            if not isinstance(phases, list) or not any(isinstance(p, dict) and p.get("name") == "setup"
                    and p.get("status") in {"ok", "completed"} for p in phases):
                raise StudyError("Setup has not recorded successful completion. Simulation is waiting.", code="setup.review.unavailable")
            if any(isinstance(p, dict) and p.get("name") == "simulation" for p in phases):
                raise StudyError("This study already records a simulation. It cannot be accepted again.", code="setup.review.unavailable")
            # Even a failed first attempt leaves records: never restart it
            # under an old acceptance or silently overwrite its output.
            simulation = root / "simulation"
            from fastmdxplora.record_access import is_dir, iterdir
            if is_dir(simulation) and any(iterdir(simulation)):
                raise StudyError("Simulation output is already present. This review cannot start it again.", code="setup.review.unavailable")
            audit, _ = _audit(root)
            if audit.get("status") != "complete":
                raise StudyError("The preparation audit is incomplete. Inspect its warnings before preparing another study.", code="setup.review.unavailable")
            final = [e for e in audit["events"] if isinstance(e, dict) and e.get("operation") == "prepared_system"]
            if len(final) != 1:
                raise StudyError("The audit does not identify one final prepared system.", code="setup.review.unavailable")
            for source_id in audit["sources"]:
                snapshot, reason = _snapshot(root, audit["sources"], source_id)
                if snapshot is None:
                    raise StudyError(reason, code="setup.review.unavailable")
                hashes[snapshot["source"]] = _hash(root / snapshot["source"])
            source = audit["sources"].get(final[0].get("after"), {})
            if source.get("sha256") != hashes["setup/topology.pdb"]:
                raise StudyError("The prepared topology no longer matches the final audited structure.", code="setup.review.changed")
            import yaml
            plan = yaml.safe_load(_read(root / PLAN, 2_000_000))
            simulation_settings = {**(plan.get("simulation") or {}), **(plan["systems"][0].get("simulation") or {})}
            if scope.denied:
                raise OutsideWorkspace("A preparation source was not admitted.")
        revision = hashlib.sha256(json.dumps(hashes, sort_keys=True).encode()).hexdigest()
        if record.get("setup_revision") and record["setup_revision"] != revision:
            raise StudyError("Prepared files changed after this review was recorded. Prepare a new audited study before accepting a different system.",
                             code="setup.review.changed")
        if not record.get("setup_revision"):
            record["setup_revision"] = revision
            _write(root, record)
        return {"requested": True, "can_accept": True, "status": "waiting", "revision": revision,
                **({"sources": hashes} if sources else {}),
                "simulation": simulation_settings,
                "reason": "Setup is complete. Simulation is waiting for your acceptance."}
    except StudyError as exc:
        return {"requested": True, "can_accept": False, "status": "waiting",
                "reason": str(exc), "refusal": exc.refusal.as_dict()}
    except (OSError, ValueError, TypeError, KeyError, UnicodeError, OutsideWorkspace):
        return {"requested": True, "can_accept": False, "status": "waiting",
                "reason": "The setup evidence or simulation plan is missing, changed, incomplete or outside this study. Refresh after setup finishes; inspect the preparation record if this persists."}


def _freeze(root: Path, hashes: dict[str, str]) -> Path:
    """Copy the checked inputs so later edits cannot change the child inputs."""
    import yaml

    token = uuid.uuid4().hex
    kept = root / ("accepted_setup-" + token)
    launch_plan = root / ("accepted_simulation-" + token + ".yml")
    kept.mkdir()
    try:
        with confined_records(root, lambda path: _public(root, path)):
            for name, expected in hashes.items():
                source = root / name
                data = _read(source, MAX_CORE).encode("utf-8")
                if hashlib.sha256(data).hexdigest() != expected:
                    raise StudyError("The preparation changed while retaining the accepted inputs.", code="setup.review.changed")
                if name == PLAN:
                    plan = yaml.safe_load(data)
                if name.startswith("setup/"):
                    target = kept / name
                    target.parent.mkdir(parents=True, exist_ok=True)
                    target.write_bytes(data)
        # Only the source location changes. The scientific settings stay
        # exactly as retained, including per-system overrides.
        location = str(kept / "setup")
        plan.setdefault("simulation", {})["setup_from"] = location
        plan["systems"][0].setdefault("simulation", {})["setup_from"] = location
        # Keep the plan beside the original so relative input paths retain
        # their meaning to the existing Config loader.
        with launch_plan.open("x", encoding="utf-8") as stream:
            stream.write(yaml.safe_dump(plan, sort_keys=False))
        launch_plan.chmod(0o444)
        for path in kept.rglob("*"):
            if path.is_file():
                path.chmod(0o444)
        return launch_plan
    except Exception:
        shutil.rmtree(kept)
        launch_plan.unlink(missing_ok=True)
        raise


def accept(root: Path, revision: str, spawn: Any, *, preflight: Any = None) -> dict[str, Any]:
    """Called with the runtime and workspace launch locks already held."""
    current = review_state(root, sources=True)
    if not revision or not current.get("can_accept") or current.get("revision") != revision:
        return {"ok": False, "error": (current.get("reason") or "This study has no pending preparation review.") if not current.get("can_accept")
                else "The preparation changed since this review. Refresh and inspect it before accepting."}
    if preflight is not None:
        refused = preflight()
        if refused is not None:
            return refused
        # Recheck after discovery, including evidence changed while probing.
        return accept(root, revision, spawn)
    try:
        plan = _freeze(root, current["sources"])
    except StudyError as exc:
        return {"ok": False, "error": str(exc), "code": exc.refusal.code}
    # The child reads its own retained copy, never the editable setup or
    # plan. No command, Config or path from the request is used here.
    with confined_records(root, lambda path: path in {root, root / RECORD}):
        record = json.loads(_read(root / RECORD, 64 * 1024))
    record.update(status="starting", revision=revision, launch_plan=str(plan.relative_to(root)))
    _write(root, record)
    command = [sys.executable, "-m", "fastmdxplora", "simulate", "--config", str(plan),
               "--output", str(root), "--no-defaults"]
    try:
        started = spawn(command, plan)
    except Exception:
        # A callback may raise after Popen succeeded, while registering the
        # process or preparing its response. Its outcome is unresolved;
        # never turn that uncertainty into permission to launch again.
        raise
    record["status"] = "accepted" if started.get("launched") and started.get("ok") is not False else "waiting"
    _write(root, record)
    return {**started, "ok": record["status"] == "accepted", "review": record["status"]}
