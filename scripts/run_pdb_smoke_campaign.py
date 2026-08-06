#!/usr/bin/env python
"""Run a small, resumable FastMDXplora smoke campaign across PDB inputs.

The campaign intentionally starts small. It accepts local PDB/CIF paths or
4-character PDB IDs, runs the normal FastMDXplora pipeline with the gentle
short simulation, validates expected artifacts, and writes CSV/JSON summaries
that make failures reviewable later.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import time
import traceback
from collections.abc import Iterable
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

SUMMARY_FIELDS = [
    "pdb_id",
    "input_path",
    "status",
    "failed_phase",
    "failure_category",
    "bug_classification",
    "error_type",
    "error_message",
    "setup_completed",
    "simulation_completed",
    "analysis_completed",
    "report_completed",
    "has_state_minimized",
    "has_trajectory",
    "runtime_seconds",
    "output_dir",
]

PHASES = ("setup", "simulation", "analysis", "report")


@dataclass
class ValidationResult:
    """Validation findings for one campaign item."""

    issues: list[str] = field(default_factory=list)
    has_state_minimized: bool = False
    has_trajectory: bool = False
    particle_count_consistent: bool | None = None

    @property
    def ok(self) -> bool:
        return not self.issues


def parse_inputs(args: argparse.Namespace) -> list[str]:
    """Return inputs from positional args and input-list files, preserving order."""
    values: list[str] = []
    for list_file in args.input_list or []:
        values.extend(_read_input_list(Path(list_file)))
    values.extend(args.inputs or [])
    return _dedupe_preserve_order(values)


def _read_input_list(path: Path) -> list[str]:
    if not path.exists():
        raise FileNotFoundError(f"Input list not found: {path}")
    values: list[str] = []
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.split("#", 1)[0].strip()
        if line:
            values.append(line)
    return values


def _dedupe_preserve_order(values: Iterable[str]) -> list[str]:
    seen: set[str] = set()
    out: list[str] = []
    for value in values:
        key = value.strip()
        if key and key not in seen:
            seen.add(key)
            out.append(key)
    return out


def input_label(value: str) -> str:
    """Return a filesystem-safe label for one input."""
    p = Path(value)
    if p.exists():
        stem = p.stem or "local_pdb"
    elif len(value) == 4 and value.isalnum():
        stem = value.upper()
    else:
        stem = p.stem or value
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", stem).strip("._") or "input"


def pdb_id_for(value: str) -> str:
    """Best-effort PDB ID for summaries; local files may not have one."""
    if len(value) == 4 and value.isalnum() and not Path(value).exists():
        return value.upper()
    stem = Path(value).stem
    return stem.upper() if len(stem) == 4 and stem.isalnum() else ""


def build_simulation_options(args: argparse.Namespace) -> dict[str, Any]:
    """Short runs at the conditions a real one uses.

    There was a `gentle` preset here that also dropped the temperature to
    100 K. A campaign wants to know whether these structures prepare and
    simulate, and a run at 100 K cannot answer that: water is ice there, and
    surviving it says nothing about 300 K.
    """
    return {
        "nvt_steps": args.nvt_steps,
        "npt_steps": args.npt_steps,
        "production_steps": args.production_steps,
        "trajectory_interval_steps": args.trajectory_interval_steps,
        "platform": args.platform,
    }


def run_one(input_value: str, args: argparse.Namespace) -> dict[str, Any]:
    """Run one FastMDXplora campaign item and return a summary row."""
    from fastmdxplora import FastMDXplora

    label = input_label(input_value)
    output_dir = Path(args.output_root) / label
    output_dir.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    row = _base_row(input_value, output_dir)
    too_large = _local_input_too_large(input_value, args.max_input_mb)
    if too_large is not None:
        row.update(
            {
                "status": "skipped",
                "failure_category": "too large for smoke settings",
                "error_type": "InputSizeLimit",
                "error_message": too_large,
            }
        )
        row["runtime_seconds"] = f"{time.monotonic() - started:.3f}"
        row["bug_classification"] = classify_bug_likelihood(row)
        return row

    try:
        fmdx = FastMDXplora(
            system=input_value,
            output_dir=output_dir,
            verbose=bool(args.verbose),
        )
        if args.max_setup_atoms and args.max_setup_atoms > 0:
            setup_results = fmdx.explore(include=["setup"])
            setup_result = setup_results[0] if setup_results else None
            setup_summary = _summarize_run_result(setup_result, validate_output_dir(output_dir))
            if setup_summary["status"] not in {"ok", "validation_failed"}:
                row.update(setup_summary)
                row["runtime_seconds"] = f"{time.monotonic() - started:.3f}"
                if not row["bug_classification"]:
                    row["bug_classification"] = classify_bug_likelihood(row)
                return row
            atom_count = _topology_atom_count(output_dir / "setup" / "topology.pdb")
            if atom_count is not None and atom_count > args.max_setup_atoms:
                row.update(setup_summary)
                row.update(
                    {
                        "status": "skipped",
                        "failed_phase": "simulation",
                        "failure_category": "too large for smoke settings",
                        "bug_classification": "expected limitation/input issue",
                        "error_type": "PreparedSystemSizeLimit",
                        "error_message": (
                            f"Prepared system has {atom_count} atoms, above "
                            f"--max-setup-atoms={args.max_setup_atoms}; "
                            "skip simulation in this smoke campaign."
                        ),
                    }
                )
                row["runtime_seconds"] = f"{time.monotonic() - started:.3f}"
                return row
            phases = ["simulation", "analysis"]
            if not args.no_report:
                phases.append("report")
            results = fmdx.explore(
                include=phases,
                options={"simulation": build_simulation_options(args)},
                report=not args.no_report,
            )
        else:
            results = fmdx.explore(
                options={"simulation": build_simulation_options(args)},
                report=not args.no_report,
            )
        run_result = results[0] if results else None
        validation = validate_output_dir(output_dir)
        row.update(_summarize_run_result(run_result, validation))
    except Exception as exc:  # noqa: BLE001 -- campaign must continue and summarize
        (output_dir / "campaign_exception.txt").write_text(
            traceback.format_exc(), encoding="utf-8"
        )
        row.update(
            {
                "status": "failed",
                "failed_phase": infer_phase_from_message(str(exc)),
                "failure_category": classify_failure(str(exc), exc_type=type(exc).__name__),
                "error_type": type(exc).__name__,
                "error_message": str(exc),
            }
        )

    row["runtime_seconds"] = f"{time.monotonic() - started:.3f}"
    if not row["bug_classification"]:
        row["bug_classification"] = classify_bug_likelihood(row)
    return row


def _base_row(input_value: str, output_dir: Path) -> dict[str, Any]:
    path = Path(input_value)
    return {
        "pdb_id": pdb_id_for(input_value),
        "input_path": str(path.resolve()) if path.exists() else input_value,
        "status": "unknown",
        "failed_phase": "",
        "failure_category": "unknown",
        "bug_classification": "",
        "error_type": "",
        "error_message": "",
        "setup_completed": False,
        "simulation_completed": False,
        "analysis_completed": False,
        "report_completed": False,
        "has_state_minimized": False,
        "has_trajectory": False,
        "runtime_seconds": "0.000",
        "output_dir": str(output_dir),
    }


def _local_input_too_large(input_value: str, max_input_mb: float | None) -> str | None:
    if max_input_mb is None or max_input_mb <= 0:
        return None
    path = Path(input_value)
    if not path.exists() or not path.is_file():
        return None
    size_mb = path.stat().st_size / (1024 * 1024)
    if size_mb <= max_input_mb:
        return None
    return (
        f"Local input is {size_mb:.1f} MB, above --max-input-mb={max_input_mb:g}; "
        "skip or run separately with a larger smoke-test budget."
    )


def _topology_atom_count(topology: Path) -> int | None:
    try:
        return sum(
            1
            for line in topology.read_text(encoding="utf-8").splitlines()
            if line.startswith(("ATOM  ", "HETATM"))
        )
    except OSError:
        return None


def _summarize_run_result(run_result: Any, validation: ValidationResult) -> dict[str, Any]:
    phases = {p.name: p for p in getattr(run_result, "phases", [])} if run_result else {}
    errored = next((p for p in phases.values() if p.status == "error"), None)
    messages = collect_messages(Path(getattr(run_result, "output_dir", "")) if run_result else None)
    status = "ok"
    failed_phase = ""
    error_type = ""
    error_message = ""

    message_text = "; ".join(messages)
    message_category = classify_failure(message_text, exc_type="RecordedPhaseNote")

    if errored is not None:
        status = "failed"
        failed_phase = errored.name
        error_message = errored.message
        error_type = infer_error_type(error_message)
    elif validation.issues and _is_expected_limitation_category(message_category):
        status = "expected_limitation"
        failed_phase = infer_phase_from_message(message_text)
        error_message = "; ".join([message_text, *validation.issues])
        error_type = "RecordedPhaseNote"
    elif validation.issues:
        status = "validation_failed"
        failed_phase = infer_phase_from_validation(validation.issues)
        error_message = "; ".join(validation.issues)
        error_type = "ValidationError"
    elif messages:
        status = "expected_limitation"
        failed_phase = infer_phase_from_message(message_text)
        error_message = message_text
        error_type = "RecordedPhaseNote"

    category = classify_failure(error_message, exc_type=error_type)
    row = {
        "status": status,
        "failed_phase": failed_phase,
        "failure_category": category if status != "ok" else "",
        "error_type": error_type,
        "error_message": error_message,
        "setup_completed": _phase_ok(phases, "setup"),
        "simulation_completed": _phase_ok(phases, "simulation"),
        "analysis_completed": _phase_ok(phases, "analysis"),
        "report_completed": _phase_ok(phases, "report"),
        "has_state_minimized": validation.has_state_minimized,
        "has_trajectory": validation.has_trajectory,
    }
    row["bug_classification"] = classify_bug_likelihood(row)
    return row


def _phase_ok(phases: dict[str, Any], phase: str) -> bool:
    result = phases.get(phase)
    return bool(result and result.status == "ok")


def validate_output_dir(output_dir: Path) -> ValidationResult:
    """Validate artifacts that should exist after reported phase success."""
    result = ValidationResult()
    manifest = _read_json(output_dir / "manifest.json")
    phase_status = {
        p.get("name"): p.get("status")
        for p in manifest.get("phases", [])
        if isinstance(p, dict)
    }

    setup_dir = output_dir / "setup"
    sim_dir = output_dir / "simulation"
    if phase_status.get("setup") == "ok":
        _require_files(
            result,
            setup_dir,
            ["setup_parameters.json", "system.xml", "state.xml", "topology.pdb"],
            "setup",
        )
    if phase_status.get("simulation") == "ok":
        _require_files(
            result,
            sim_dir,
            ["simulation_parameters.json", "state_minimized.xml", "production.dcd"],
            "simulation",
        )
    if phase_status.get("analysis") == "ok":
        _require_files(result, output_dir / "analysis", ["analysis_manifest.json"], "analysis")
    if phase_status.get("report") == "ok":
        _require_files(result, output_dir / "report", ["report.md"], "report")

    result.has_state_minimized = (sim_dir / "state_minimized.xml").exists()
    result.has_trajectory = (sim_dir / "production.dcd").exists()
    result.issues.extend(check_pdb_coordinates(output_dir))
    result.issues.extend(
        check_trajectory_coordinates(sim_dir / "production.dcd", sim_dir / "topology.pdb")
    )
    particle_check = check_particle_counts(output_dir)
    if particle_check is not None:
        result.particle_count_consistent = particle_check[0]
        if not particle_check[0]:
            result.issues.append(particle_check[1])
    return result


def _require_files(result: ValidationResult, base: Path, names: list[str], phase: str) -> None:
    for name in names:
        if not (base / name).exists():
            result.issues.append(f"{phase}: missing expected artifact {base / name}")


def check_pdb_coordinates(output_dir: Path) -> list[str]:
    """Check saved PDB coordinate columns for NaN/Inf."""
    issues: list[str] = []
    for pdb in output_dir.rglob("*.pdb"):
        try:
            for line_no, line in enumerate(pdb.read_text(encoding="utf-8").splitlines(), 1):
                if not line.startswith(("ATOM  ", "HETATM")):
                    continue
                xyz = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
                if not all(math.isfinite(v) for v in xyz):
                    issues.append(f"{pdb}: non-finite PDB coordinate on line {line_no}")
                    break
        except (OSError, ValueError) as exc:
            issues.append(f"{pdb}: could not validate PDB coordinates ({exc})")
    return issues


def check_trajectory_coordinates(trajectory: Path, topology: Path) -> list[str]:
    """Best-effort trajectory NaN/Inf check using MDTraj when available."""
    if not trajectory.exists() or not topology.exists():
        return []
    try:
        import mdtraj as md
        import numpy as np
    except ImportError:
        return []
    try:
        traj = md.load(str(trajectory), top=str(topology))
    except Exception as exc:  # noqa: BLE001 -- bad trajectory is validation signal
        return [f"{trajectory}: could not load trajectory for validation ({exc})"]
    if not np.isfinite(traj.xyz).all():
        return [f"{trajectory}: non-finite trajectory coordinate detected"]
    return []


def check_particle_counts(output_dir: Path) -> tuple[bool, str] | None:
    """Compare OpenMM System particle count with topology PDB atom count when possible."""
    system_xml = output_dir / "setup" / "system.xml"
    topology = output_dir / "setup" / "topology.pdb"
    if not system_xml.exists() or not topology.exists():
        return None
    try:
        n_particles = _openmm_system_particle_count(system_xml)
        if n_particles is None:
            text = system_xml.read_text(encoding="utf-8")
            n_particles = len(re.findall(r"<Particle\s+mass=", text))
        n_atoms = sum(
            1
            for line in topology.read_text(encoding="utf-8").splitlines()
            if line.startswith(("ATOM  ", "HETATM"))
        )
    except OSError:
        return None
    if n_particles and n_atoms and n_particles != n_atoms:
        return False, f"particle count mismatch: system.xml={n_particles}, topology.pdb={n_atoms}"
    return True, ""


def _openmm_system_particle_count(system_xml: Path) -> int | None:
    try:
        from openmm import XmlSerializer
    except ImportError:
        return None
    try:
        system = XmlSerializer.deserialize(system_xml.read_text(encoding="utf-8"))
        return int(system.getNumParticles())
    except Exception:  # noqa: BLE001 -- fall back to text validation
        return None


def collect_messages(output_dir: Path | None) -> list[str]:
    """Collect non-empty phase manifest notes that indicate deferred/skipped work."""
    if output_dir is None or not output_dir.exists():
        return []
    messages: list[str] = []
    for manifest in (
        output_dir / "setup" / "setup_parameters.json",
        output_dir / "simulation" / "simulation_parameters.json",
        output_dir / "analysis" / "analysis_manifest.json",
    ):
        data = _read_json(manifest)
        for note in data.get("notes", []) or []:
            if note:
                messages.append(str(note))
        if data.get("status") == "deferred" and data.get("note"):
            messages.append(str(data["note"]))
    return messages


def _read_json(path: Path) -> dict[str, Any]:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}


def classify_failure(message: str, *, exc_type: str = "") -> str:
    """Classify an error or validation message into a review category."""
    text = f"{exc_type} {message}".lower()
    if (
        "no such file" in text
        or "not found" in text and "setup outputs" not in text
        or "could not classify system input" in text and ".pdb" in text
    ):
        return "missing input file"
    if any(s in text for s in ("temporary failure in name resolution", "dns", "urlopen")):
        return "DNS/download failure"
    if "openmm unavailable" in text or "pdbfixer unavailable" in text:
        return "missing dependency"
    if any(s in text for s in ("template", "no template", "residue", "nonstandard")):
        return "unsupported residue/template failure"
    if "ligand" in text or "openff" in text or "rdkit" in text:
        return "ligand unsupported"
    if any(s in text for s in ("metal", "zn", "mg", "fe", "heme")):
        return "metal unsupported"
    if any(s in text for s in ("missing atom", "missing residue", "cannot be repaired")):
        return "missing atoms/residues issue"
    if any(s in text for s in ("geometry", "clash", "nan", "particle coordinate is nan")):
        return "bad geometry/clash" if "nan" not in text else "OpenMM NaN"
    if any(s in text for s in ("box", "solvat", "periodic box")):
        return "solvation/box issue"
    if "too large" in text or "above --max-input-mb" in text:
        return "too large for smoke settings"
    if "analysis" in text:
        return "analysis failure"
    if "report" in text or "pptx" in text:
        return "report generation failure"
    if exc_type and exc_type not in {"RecordedPhaseNote", "ValidationError"}:
        return "code exception/bug"
    return "unknown"


def classify_bug_likelihood(row: dict[str, Any]) -> str:
    """Apply campaign policy for code bugs vs expected input/environment limitations."""
    category = row.get("failure_category", "")
    status = row.get("status", "")
    message = str(row.get("error_message", "")).lower()
    if status == "ok":
        return "none"
    if status == "validation_failed":
        return "likely code bug"
    if _is_expected_limitation_category(category):
        return "expected limitation/input issue"
    if "internal" in message or category == "code exception/bug":
        return "likely code bug"
    if category in {"analysis failure", "report generation failure", "OpenMM NaN"}:
        return "likely code bug"
    return "needs review"


def _is_expected_limitation_category(category: str) -> bool:
    expected = {
        "DNS/download failure",
        "unsupported residue/template failure",
        "ligand unsupported",
        "metal unsupported",
        "missing atoms/residues issue",
        "solvation/box issue",
        "missing dependency",
        "missing input file",
        "too large for smoke settings",
    }
    return category in expected


def infer_error_type(message: str) -> str:
    match = re.match(r"([A-Za-z_][A-Za-z0-9_]*):", message.strip())
    return match.group(1) if match else ""


def infer_phase_from_message(message: str) -> str:
    text = message.lower()
    for phase in PHASES:
        if phase in text:
            return phase
    if "pdbfixer" in text or "template" in text or "download" in text:
        return "setup"
    if "openmm" in text or "nan" in text or "trajectory" in text:
        return "simulation"
    return ""


def infer_phase_from_validation(issues: list[str]) -> str:
    joined = " ".join(issues).lower()
    for phase in PHASES:
        if f"{phase}:" in joined:
            return phase
    return infer_phase_from_message(joined)


def write_summaries(rows: list[dict[str, Any]], output_root: Path) -> None:
    output_root.mkdir(parents=True, exist_ok=True)
    csv_path = output_root / "campaign_summary.csv"
    json_path = output_root / "campaign_summary.json"
    summary_csv_path = output_root / "summary.csv"
    summary_json_path = output_root / "summary.json"
    with csv_path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=SUMMARY_FIELDS)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in SUMMARY_FIELDS})
    csv_text = csv_path.read_text(encoding="utf-8")
    json_text = json.dumps(rows, indent=2, default=str)
    summary_csv_path.write_text(csv_text, encoding="utf-8")
    json_path.write_text(json_text, encoding="utf-8")
    summary_json_path.write_text(json_text, encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run a staged FastMDXplora PDB smoke campaign.",
    )
    parser.add_argument("inputs", nargs="*", help="PDB IDs or local .pdb/.cif files.")
    parser.add_argument(
        "--input-list",
        action="append",
        help="Text file of PDB IDs or local paths, one per line. # comments are allowed.",
    )
    parser.add_argument("--output-root", required=True, help="Directory for campaign outputs.")
    parser.add_argument("--nvt-steps", type=int, default=1000)
    parser.add_argument("--npt-steps", type=int, default=0)
    parser.add_argument("--production-steps", type=int, default=1000)
    parser.add_argument("--trajectory-interval-steps", type=int, default=100)
    parser.add_argument("--platform", default="auto")
    parser.add_argument(
        "--max-input-mb",
        type=float,
        default=10.0,
        help="Skip local input files larger than this size; <=0 disables the guard.",
    )
    parser.add_argument(
        "--max-setup-atoms",
        type=int,
        default=0,
        help=(
            "After setup, skip simulation for prepared systems above this atom count; "
            "0 disables the guard."
        ),
    )
    parser.add_argument("--no-report", action="store_true", help="Skip report generation.")
    parser.add_argument(
        "--continue-on-error",
        action="store_true",
        help="Accepted for explicitness; campaigns continue on errors by default.",
    )
    parser.add_argument(
        "--stop-on-error", action="store_true", help="Stop after the first failure."
    )
    parser.add_argument("--verbose", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    inputs = parse_inputs(args)
    if not inputs:
        parser.error("provide at least one PDB ID/local file or --input-list")

    rows: list[dict[str, Any]] = []
    output_root = Path(args.output_root)
    for value in inputs:
        row = run_one(value, args)
        rows.append(row)
        write_summaries(rows, output_root)
        print(
            f"{value}: {row['status']} "
            f"{row['failure_category'] or 'ok'} -> {row['output_dir']}",
            flush=True,
        )
        if args.stop_on_error and row["status"] != "ok":
            break
    return 0 if all(row["status"] == "ok" for row in rows) else 1


if __name__ == "__main__":
    raise SystemExit(main())
