"""Confined inspection of retained force evidence, without evaluating forces."""
from __future__ import annotations

from copy import deepcopy
import hashlib
import json
import math
import re
from pathlib import Path
from typing import Any

from fastmdxplora.gui.preparation_review import _public, _read, pdb_atoms
from fastmdxplora.record_access import confined_records, is_file
from fastmdxplora.refusals import OutsideWorkspace

RECORD = "simulation/instability.json"
MAX_RECORD = 100_000
MAX_SNAPSHOT = 32 * 1024 * 1024
UNAVAILABLE = "No retained finite force evidence is available for this run. Older runs did not record it; forces will not be recomputed."


def evidence_runs(root: Path) -> list[tuple[Path, str]]:
    """Only declared runs; selecting evidence never reads every run's status."""
    if not is_file(root / "batch_manifest.json"):
        return [(root, "This run")]
    manifest = json.loads(_read(root / "batch_manifest.json", MAX_RECORD))
    if not isinstance(manifest, dict):
        return []
    entries = manifest.get("planned") or manifest.get("runs") or []
    if not isinstance(entries, list) or len(entries) > 128:
        return []
    result = []
    seen = set()
    for entry in entries:
        if not isinstance(entry, dict):
            return []
        name = entry.get("run_id")
        if not isinstance(name, str) or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,127}", name) or name in seen:
            return []
        seen.add(name)
        result.append((root / "runs" / name, name))
    return result


def _finite(value: Any) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)


def _evidence(root: Path) -> dict[str, Any]:
    raw = _read(root / RECORD, MAX_RECORD)
    record = json.loads(raw)
    if not isinstance(record, dict) or type(record.get("version")) is not int or record["version"] != 1:
        return {"available": False, "reason": "The force evidence format is unsupported."}
    if record.get("available") is not True:
        return {"available": False, "reason": str(record.get("reason") or UNAVAILABLE)[:1000]}
    name, digest = record.get("snapshot"), record.get("snapshot_sha256")
    if not isinstance(digest, str) or not re.fullmatch(r"[0-9a-f]{64}", digest) or name != f"instability-snapshot-{digest[:20]}.pdb":
        return {"available": False, "reason": "The force evidence does not name a supported snapshot."}
    snapshot = _read(root / "simulation" / name, MAX_SNAPSHOT)
    if hashlib.sha256(snapshot.encode("utf-8")).hexdigest() != digest:
        return {"available": False, "reason": "The retained snapshot changed. No force ranking will be mapped to it."}
    atoms, reason = pdb_atoms(snapshot)
    if reason or not atoms or len(atoms) != record.get("atom_count"):
        return {"available": False, "reason": reason or "The snapshot's full-system atom count does not match."}
    ranking = record.get("ranking")
    if not isinstance(ranking, list) or not 1 <= len(ranking) <= 12:
        return {"available": False, "reason": "No bounded force ranking can be read."}
    indices = set()
    previous = math.inf
    for row in ranking:
        if not isinstance(row, dict):
            return {"available": False, "reason": "The force ranking contains an unreadable atom."}
        index = row.get("index")
        if not isinstance(index, int) or isinstance(index, bool) or not 0 <= index < len(atoms) or index in indices:
            return {"available": False, "reason": "The full-system atom correspondence is invalid."}
        indices.add(index)
        atom = atoms[index]
        xyz, vector, norm = row.get("position_nm"), row.get("vector_kj_mol_nm"), row.get("force_kj_mol_nm")
        if not all(isinstance(values, list) and len(values) == 3 and all(_finite(v) for v in values) for values in (xyz, vector)):
            return {"available": False, "reason": "A force or coordinate is not a finite vector."}
        if not _finite(norm) or norm < 0 or norm > previous or not math.isclose(math.hypot(*vector), norm, rel_tol=1e-9, abs_tol=1e-9):
            return {"available": False, "reason": "The force ranking has inconsistent magnitudes."}
        previous = norm
        if row.get("identity") != atom["identity"] or any(abs(xyz[i] * 10 - atom["xyz"][i]) > .0006 for i in range(3)):
            return {"available": False, "reason": "A ranked atom does not match the saved identity and coordinates."}
    for key in ("step", "failure_reported_step", "interval_steps"):
        if not isinstance(record.get(key), int) or isinstance(record[key], bool) or record[key] < 0:
            return {"available": False, "reason": "The observation's step or interval is invalid."}
    if record["failure_reported_step"] < record["step"] or not _finite(record.get("timestep_fs")) or record["timestep_fs"] <= 0:
        return {"available": False, "reason": "The observation's timestep or failure counter is inconsistent."}
    if record["interval_steps"] < 1 or not _finite(record.get("time_ps")) or record["time_ps"] < 0 or record.get("force_units") != "kJ mol^-1 nm^-1" or record.get("position_units") != "nm":
        return {"available": False, "reason": "The observation's time, units or interval is unsupported."}
    if record.get("stage") != record.get("failure_stage") or record.get("stage") not in {"minimization", "nvt", "npt", "production"}:
        return {"available": False, "reason": "The finite observation does not belong to the failing stage."}
    version = hashlib.sha256(raw.encode("utf-8")).hexdigest()
    # Send only declared fields; a hand-edited file cannot inject arbitrary
    # values into the page's evidence or proposal.
    return {"available": True, "version": version, "pdb": snapshot,
        "snapshot_sha256": digest, "source": "simulation/" + name,
        "ranking": [{key: row[key] for key in ("index", "identity", "position_nm", "force_kj_mol_nm", "vector_kj_mol_nm")} for row in ranking],
        "atom_count": len(atoms), "stage": record["stage"], "step": record["step"],
        "time_ps": record["time_ps"], "failure_reported_step": record["failure_reported_step"],
        "interval_steps": record["interval_steps"], "force_units": record["force_units"],
        "reported_step_gap": record["failure_reported_step"] - record["step"],
        "notice": "Last recorded finite observation, not necessarily the step before failure. Large forces do not establish the cause.",
        "failure_step_note": "The reported step counter may precede the failed integrator step; the exact failing step is unknown.",
        "force_scope": "Active integration force groups, including their restraints and bias.",
        "system_sha256": str(record.get("system_sha256") or "")[:64],
        "config_sha256": str(record.get("config_sha256") or "")[:64],
        "timestep_fs": record.get("timestep_fs"),
        "failure": str(record.get("failure") or "")[:4000]}


def instability_payload(root: str | Path, *, boundary: Path | None = None) -> dict[str, Any]:
    root = Path(root).absolute()
    admitted = Path(boundary).absolute() if boundary is not None else root.resolve()
    try:
        with confined_records(admitted, lambda path: _public(admitted, path)) as scope:
            answer = _evidence(root)
        if scope.denied:
            return {"available": False, "reason": "A force evidence source is outside the study's public files."}
        return answer
    except (OSError, ValueError, TypeError, KeyError, OverflowError, UnicodeError, OutsideWorkspace):
        return {"available": False, "reason": UNAVAILABLE}


def timestep_proposal(root: str | Path, version: str, *, boundary: Path | None = None) -> dict[str, Any]:
    """A validated new Config for a controlled integration test, never a cure.

    Only an ordinary fixed-step single run is supported. Step-based durations
    and output intervals are doubled to keep their physical times unchanged.
    The original Config and its provenance are preserved in full.
    """
    root = Path(root).absolute()
    admitted = Path(boundary).absolute() if boundary is not None else root.resolve()
    evidence = instability_payload(root, boundary=admitted)
    if not evidence.get("available") or evidence.get("version") != version:
        return {"ok": False, "reason": "Read the current finite force evidence before requesting a proposal."}
    if evidence["stage"] == "minimization":
        return {"ok": False, "reason": "A dynamics timestep cannot fix a minimization failure. Inspect the ranked atoms and preparation."}
    try:
        import yaml
        from fastmdxplora.gui.config_yaml_edit import _StrictLoader, _bounded, _guard_values
        from fastmdxplora.config.loader import canonical_phase_keys, validate_config

        with confined_records(admitted, lambda path: _public(admitted, path)) as scope:
            raw = _read(root / "resolved_config.yml", 1_000_000)
        if scope.denied:
            return {"ok": False, "reason": "The run Config is outside this study's public files."}
        if not evidence["config_sha256"] or hashlib.sha256(raw.encode()).hexdigest() != evidence["config_sha256"]:
            return {"ok": False, "reason": "The saved Config is missing or changed since the force observation. No proposal was inferred from it."}
        _bounded(raw)
        original = yaml.load(raw, Loader=_StrictLoader)
        _guard_values(original)
        validate_config(deepcopy(original), require_systems=True)
        config = canonical_phase_keys(deepcopy(original))
        sim = config.get("simulation") or {}
        systems = config.get("systems") or []
        if not isinstance(sim, dict) or len(systems) != 1 or config.get("sweep") or any(
                sim.get(key) for key in ("umbrella", "steered", "metadynamics", "plumed", "stop_when", "resume_from", "prepared_from", "setup_from")) or any(
                isinstance(item, dict) and item.get("simulation") for item in systems):
            return {"ok": False, "reason": "This proposal supports one unbiased run without overrides, adaptive stopping or a continuation. Keep its specialized plan and inspect preparation."}
        dt = sim.get("timestep_fs")
        if not _finite(dt) or dt <= 0 or dt / 2 <= 0 or dt != evidence.get("timestep_fs") or str(sim.get("integrator", "langevin_middle")).lower().startswith("variable"):
            return {"ok": False, "reason": "The run Config does not declare a supported fixed timestep."}
        changes = [{"setting": "simulation.timestep_fs", "before": dt, "after": dt / 2}]
        sim["timestep_fs"] = dt / 2
        from fastmdxplora.config.schema import PHASE_SCHEMAS

        defaults = {field.name: field.default for field in PHASE_SCHEMAS["simulation"].fields}
        for key in ("state_interval_steps", "checkpoint_interval_steps", "telemetry_interval", "barostat_frequency"):
            if sim.get(key) is None:
                sim[key] = defaults[key]
        for key in ("nvt_steps", "npt_steps", "production_steps", "trajectory_interval_steps", "state_interval_steps", "checkpoint_interval_steps", "telemetry_interval", "barostat_frequency"):
            if key not in sim or sim[key] is None:
                continue
            value = sim[key]
            if not isinstance(value, int) or isinstance(value, bool) or value < 0:
                return {"ok": False, "reason": "A step-based duration or interval cannot be preserved safely."}
            sim[key] = value * 2
            if value:
                changes.append({"setting": "simulation." + key, "before": value, "after": value * 2})
        config["simulation"] = sim
        old_output = config.get("output")
        config["output"] = str(old_output).rstrip("/") + "_timestep_test" if old_output else "instability_timestep_test"
        changes.append({"setting": "output", "before": old_output, "after": config["output"]})
        decisions = deepcopy(config.get("decisions") or {})
        for change in changes:
            name = change["setting"]
            earlier = decisions.get(name)
            decisions[name] = {"source": "FastMDXplora instability proposal (unconfirmed)",
                "why": "Proposed controlled timestep test from finite force evidence " + version + "; not a causal diagnosis or an accepted decision.",
                "alternatives": [{"value": change["before"], "earlier_decision": earlier}]}
        config["decisions"] = decisions
        _guard_values(config)
        validate_config(deepcopy(config), require_systems=True)
        return {"ok": True, "version": version, "config_sha256": hashlib.sha256(raw.encode()).hexdigest(),
            "changes": changes, "yaml": yaml.safe_dump(config, sort_keys=False),
            "reason": "Proposal: test half the fixed timestep while preserving stated durations and step-based output times. The force ranking does not prove an integration error or that this change will fix it. Inspect the ranked atoms and preparation first. A new study is required; no run has started."}
    except Exception:
        return {"ok": False, "reason": "The saved run Config could not produce a valid proposal. No settings were changed."}
