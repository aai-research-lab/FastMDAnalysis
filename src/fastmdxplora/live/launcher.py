"""Dashboard-first workflow launcher for FastMDXplora.

The live dashboard normally watches an existing project output directory.
This module adds a small, deliberately conservative launcher layer so the
same local server can start before a run exists, validate a configuration,
and launch the normal ``fastmdxplora.cli.main explore`` workflow in a child
process.  It does not reimplement any setup, simulation, analysis, or report
science.
"""

from __future__ import annotations

import importlib
import json
import os
import re
import subprocess
import sys
import threading
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping


_FORCEFIELDS = ("charmm36", "amber14", "amber-fb15", "amber-openff")
_PLATFORMS = ("auto", "CPU", "CUDA", "OpenCL", "HIP")
_PRECISIONS = ("single", "mixed", "double")
_INTEGRATORS = (
    "langevin_middle",
    "langevin",
    "brownian",
    "verlet",
    "variable_langevin",
    "variable_verlet",
)
_ANALYSES = (
    "rmsd",
    "rmsf",
    "rg",
    "hbonds",
    "sasa",
    "ss",
    "qvalue",
    "cluster",
    "dimred",
    "dihedrals",
)


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _slug(value: str, fallback: str = "fastmdxplora_run") -> str:
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", str(value).strip()).strip("._-")
    return cleaned[:96] or fallback


def _number(
    payload: Mapping[str, Any],
    key: str,
    *,
    default: float,
    minimum: float,
    maximum: float,
    integer: bool = False,
) -> int | float:
    raw = payload.get(key, default)
    try:
        value = int(raw) if integer else float(raw)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{key} must be a number") from exc
    if value < minimum or value > maximum:
        raise ValueError(f"{key} must be between {minimum:g} and {maximum:g}")
    return value


def launcher_defaults() -> dict[str, Any]:
    """Return UI defaults from the same Python sources used by the CLI."""
    from fastmdxplora.setup.pipeline import DEFAULTS as setup_defaults
    from fastmdxplora.simulation.pipeline import DEFAULTS as sim_defaults
    from fastmdxplora.simulation.runner import (
        DEFAULT_NPT_STEPS,
        DEFAULT_NVT_STEPS,
        DEFAULT_PRODUCTION_STEPS,
        DEFAULT_TRAJECTORY_INTERVAL_STEPS,
    )

    return {
        "system": "1L2Y",
        "run_name": "",
        "setup": {
            "ph": float(setup_defaults.get("ph", 7.0)),
            "forcefield": str(setup_defaults.get("forcefield", "charmm36")),
            "water_model": setup_defaults.get("water_model") or "auto",
            "ion_concentration_M": float(setup_defaults.get("ion_concentration_M", 0.15)),
            "solvent_padding_nm": float(setup_defaults.get("solvent_padding_nm", 1.0)),
            "keep_heterogens": bool(setup_defaults.get("keep_heterogens", False)),
            "keep_water": bool(setup_defaults.get("keep_water", False)),
        },
        "simulation": {
            "minimize": bool(sim_defaults.get("minimize", True)),
            "nvt_steps": int(DEFAULT_NVT_STEPS),
            "npt_steps": int(DEFAULT_NPT_STEPS),
            "production_steps": int(DEFAULT_PRODUCTION_STEPS),
            "timestep_fs": float(sim_defaults.get("timestep_fs", 2.0)),
            "temperature_K": float(sim_defaults.get("temperature_K", 300.0)),
            "friction_per_ps": float(sim_defaults.get("friction_per_ps", 1.0)),
            "integrator": str(sim_defaults.get("integrator", "langevin_middle")),
            "platform": str(sim_defaults.get("platform", "auto")),
            "precision": str(sim_defaults.get("precision", "mixed")),
            "trajectory_interval_steps": int(DEFAULT_TRAJECTORY_INTERVAL_STEPS),
            "checkpoint_interval_steps": int(sim_defaults.get("checkpoint_interval_steps", 10000)),
            "telemetry_interval": int(sim_defaults.get("telemetry_interval", 1000)),
        },
        "workflow": {
            "run_analysis": True,
            "run_report": True,
            "analyses": list(_ANALYSES),
            "report_document": True,
            "report_slides": True,
            "report_bundle": True,
        },
        "choices": {
            "forcefields": list(_FORCEFIELDS),
            "platforms": list(_PLATFORMS),
            "precisions": list(_PRECISIONS),
            "integrators": list(_INTEGRATORS),
            "analyses": list(_ANALYSES),
        },
    }


def validate_launcher_payload(payload: Mapping[str, Any]) -> dict[str, Any]:
    """Validate and normalize one dashboard launch request.

    The returned mapping is safe to translate into a list-form subprocess
    command.  Validation is intentionally stricter than argparse so bad form
    values are reported next to their fields instead of failing after launch.
    """
    errors: dict[str, str] = {}
    warnings: list[str] = []

    system = str(payload.get("system") or "").strip()
    if not system:
        errors["system"] = "Enter a PDB ID, sequence, or local PDB/CIF path."
    elif len(system) > 4096:
        errors["system"] = "System input is too long."

    run_name_raw = str(payload.get("run_name") or "").strip()
    default_name = f"{_slug(system, 'system')}_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
    run_name = _slug(run_name_raw, default_name)

    setup_in = payload.get("setup") if isinstance(payload.get("setup"), Mapping) else {}
    sim_in = payload.get("simulation") if isinstance(payload.get("simulation"), Mapping) else {}
    workflow_in = payload.get("workflow") if isinstance(payload.get("workflow"), Mapping) else {}

    setup: dict[str, Any] = {}
    simulation: dict[str, Any] = {}
    workflow: dict[str, Any] = {}

    try:
        setup["ph"] = _number(setup_in, "ph", default=7.0, minimum=0.0, maximum=14.0)
    except ValueError as exc:
        errors["setup.ph"] = str(exc)
    forcefield = str(setup_in.get("forcefield") or "charmm36")
    if forcefield not in _FORCEFIELDS:
        errors["setup.forcefield"] = "Choose a supported force field."
    else:
        setup["forcefield"] = forcefield
    water_model = str(setup_in.get("water_model") or "auto").strip()
    setup["water_model"] = water_model
    try:
        setup["ion_concentration_M"] = _number(
            setup_in,
            "ion_concentration_M",
            default=0.15,
            minimum=0.0,
            maximum=5.0,
        )
    except ValueError as exc:
        errors["setup.ion_concentration_M"] = str(exc)
    try:
        setup["solvent_padding_nm"] = _number(
            setup_in,
            "solvent_padding_nm",
            default=1.0,
            minimum=0.1,
            maximum=10.0,
        )
    except ValueError as exc:
        errors["setup.solvent_padding_nm"] = str(exc)
    setup["keep_heterogens"] = bool(setup_in.get("keep_heterogens", False))
    setup["keep_water"] = bool(setup_in.get("keep_water", False))

    numeric_fields = (
        ("nvt_steps", 250000, 0, 10_000_000_000, True),
        ("npt_steps", 500000, 0, 10_000_000_000, True),
        ("production_steps", 1000000, 1, 100_000_000_000, True),
        ("timestep_fs", 2.0, 0.01, 20.0, False),
        ("temperature_K", 300.0, 1.0, 5000.0, False),
        ("friction_per_ps", 1.0, 0.0, 1000.0, False),
        ("trajectory_interval_steps", 1000, 1, 10_000_000_000, True),
        ("checkpoint_interval_steps", 10000, 0, 10_000_000_000, True),
        ("telemetry_interval", 1000, 1, 10_000_000_000, True),
    )
    for key, default, low, high, integer in numeric_fields:
        try:
            simulation[key] = _number(
                sim_in,
                key,
                default=default,
                minimum=low,
                maximum=high,
                integer=integer,
            )
        except ValueError as exc:
            errors[f"simulation.{key}"] = str(exc)

    simulation["minimize"] = bool(sim_in.get("minimize", True))
    integrator = str(sim_in.get("integrator") or "langevin_middle")
    platform = str(sim_in.get("platform") or "auto")
    precision = str(sim_in.get("precision") or "mixed")
    if integrator not in _INTEGRATORS:
        errors["simulation.integrator"] = "Choose a supported integrator."
    else:
        simulation["integrator"] = integrator
    if platform not in _PLATFORMS:
        errors["simulation.platform"] = "Choose a supported OpenMM platform."
    else:
        simulation["platform"] = platform
    if precision not in _PRECISIONS:
        errors["simulation.precision"] = "Choose a supported precision mode."
    else:
        simulation["precision"] = precision

    workflow["run_analysis"] = bool(workflow_in.get("run_analysis", True))
    workflow["run_report"] = bool(workflow_in.get("run_report", True))
    requested_analyses = workflow_in.get("analyses", list(_ANALYSES))
    if not isinstance(requested_analyses, list):
        requested_analyses = list(_ANALYSES)
    analyses = [str(item) for item in requested_analyses if str(item) in _ANALYSES]
    workflow["analyses"] = analyses
    workflow["report_document"] = bool(workflow_in.get("report_document", True))
    workflow["report_slides"] = bool(workflow_in.get("report_slides", True))
    workflow["report_bundle"] = bool(workflow_in.get("report_bundle", True))

    if workflow["run_report"] and not workflow["run_analysis"]:
        warnings.append("Reports can be generated without analysis, but they will contain fewer scientific results.")

    timestep = float(simulation.get("timestep_fs", 2.0))
    nvt_steps = int(simulation.get("nvt_steps", 0))
    npt_steps = int(simulation.get("npt_steps", 0))
    production_steps = int(simulation.get("production_steps", 0))
    trajectory_interval = int(simulation.get("trajectory_interval_steps", 1))
    telemetry_interval = int(simulation.get("telemetry_interval", 1))

    durations = {
        "nvt_ns": nvt_steps * timestep / 1_000_000.0,
        "npt_ns": npt_steps * timestep / 1_000_000.0,
        "production_ns": production_steps * timestep / 1_000_000.0,
        "total_ns": (nvt_steps + npt_steps + production_steps) * timestep / 1_000_000.0,
        "trajectory_frames": (production_steps + trajectory_interval - 1) // trajectory_interval,
        "dashboard_frames": (
            nvt_steps + npt_steps + production_steps + telemetry_interval - 1
        ) // telemetry_interval,
    }
    if durations["production_ns"] < 1.0:
        warnings.append(
            f"Production time is {durations['production_ns']:.4g} ns; this is suitable for testing, not strong scientific conclusions."
        )

    return {
        "valid": not errors,
        "errors": errors,
        "warnings": warnings,
        "config": {
            "system": system,
            "run_name": run_name,
            "setup": setup,
            "simulation": simulation,
            "workflow": workflow,
        },
        "summary": durations,
    }


def build_launcher_command(config: Mapping[str, Any], output_dir: Path) -> list[str]:
    """Translate a normalized launcher config into the canonical CLI."""
    setup = config["setup"]
    sim = config["simulation"]
    workflow = config["workflow"]
    command = [
        sys.executable,
        "-m",
        "fastmdxplora.cli.main",
        "explore",
        "--system",
        str(config["system"]),
        "--output",
        str(output_dir),
        "--setup-ph",
        str(setup["ph"]),
        "--setup-forcefield",
        str(setup["forcefield"]),
        "--setup-ion-concentration-M",
        str(setup["ion_concentration_M"]),
        "--setup-solvent-padding-nm",
        str(setup["solvent_padding_nm"]),
        "--simulate-nvt-steps",
        str(sim["nvt_steps"]),
        "--simulate-npt-steps",
        str(sim["npt_steps"]),
        "--simulate-production-steps",
        str(sim["production_steps"]),
        "--simulate-timestep-fs",
        str(sim["timestep_fs"]),
        "--simulate-temperature-K",
        str(sim["temperature_K"]),
        "--simulate-friction-per-ps",
        str(sim["friction_per_ps"]),
        "--simulate-integrator",
        str(sim["integrator"]),
        "--simulate-platform",
        str(sim["platform"]),
        "--simulate-precision",
        str(sim["precision"]),
        "--simulate-trajectory-interval-steps",
        str(sim["trajectory_interval_steps"]),
        "--simulate-checkpoint-interval-steps",
        str(sim["checkpoint_interval_steps"]),
        "--simulate-live-telemetry",
        "--simulate-telemetry-interval",
        str(sim["telemetry_interval"]),
    ]
    if setup.get("water_model") and setup["water_model"] != "auto":
        command.extend(["--setup-water-model", str(setup["water_model"])])
    if setup.get("keep_heterogens"):
        command.append("--setup-keep-heterogens")
    if setup.get("keep_water"):
        command.append("--setup-keep-water")
    if not sim.get("minimize", True):
        command.append("--simulate-no-minimize")

    phases = ["setup", "simulation"]
    if workflow.get("run_analysis"):
        phases.append("analysis")
    if workflow.get("run_report"):
        phases.append("report")
    if phases != ["setup", "simulation", "analysis", "report"]:
        command.extend(["--include", *phases])

    analyses = list(workflow.get("analyses") or [])
    if workflow.get("run_analysis") and analyses and set(analyses) != set(_ANALYSES):
        command.extend(["--analyze-analyses", *analyses])
    if workflow.get("run_report"):
        if not workflow.get("report_document", True):
            command.append("--report-no-document")
        if not workflow.get("report_slides", True):
            command.append("--report-no-slides")
        if not workflow.get("report_bundle", True):
            command.append("--report-no-bundle")
        command.extend(["--report-title", str(config.get("run_name") or "FastMDXplora Run")])
    return command


def launcher_environment_error(config: Mapping[str, Any]) -> str | None:
    """Return an actionable error when the dashboard cannot run the workflow.

    The normal CLI deliberately imports the chemistry stack lazily and lets
    phase pipelines write a warning manifest when optional packages are
    missing.  That behavior is useful for inspecting configuration on light
    installations, but it is misleading for the dashboard's *Run Simulation*
    action: the child exits successfully without producing a simulation.
    """
    required = [
        ("OpenMM", "openmm"),
        ("OpenMM application layer", "openmm.app"),
        ("PDBFixer", "pdbfixer"),
    ]
    workflow = config.get("workflow")
    if isinstance(workflow, Mapping) and workflow.get("run_analysis"):
        required.append(("MDTraj", "mdtraj"))

    missing: list[str] = []
    for label, module_name in required:
        try:
            importlib.import_module(module_name)
        except Exception:  # noqa: BLE001 - broken partial installs also cannot run
            if label.startswith("OpenMM"):
                label = "OpenMM"
            if label not in missing:
                missing.append(label)

    if not missing:
        return None

    packages = " ".join(
        name for name in ("openmm", "pdbfixer", "mdtraj")
        if name != "mdtraj" or "MDTraj" in missing
    )
    return (
        "Simulation dependencies are unavailable in the Python environment "
        f"running this dashboard: {', '.join(missing)}. Install them in that "
        "same environment, restart the dashboard, and launch again. Recommended: "
        f"conda install -c conda-forge {packages}"
    )


def _json_mapping(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {}
    return value if isinstance(value, dict) else {}


def _manifest_note(manifest: Mapping[str, Any]) -> str | None:
    notes = manifest.get("notes")
    if not isinstance(notes, list):
        return None
    for note in notes:
        text = str(note or "").strip()
        if text:
            return text
    return None


def _is_nonempty_file(path: Path) -> bool:
    try:
        return path.is_file() and path.stat().st_size > 0
    except OSError:
        return False


@dataclass
class DashboardRuntime:
    """Mutable state shared by all request-handler threads."""

    workspace_root: Path
    launch_root: Path
    active_root: Path | None = None
    process: subprocess.Popen[Any] | None = None
    process_started_at: str | None = None
    process_finished_at: str | None = None
    process_returncode: int | None = None
    completion_error: str | None = None
    log_path: Path | None = None
    command: list[str] = field(default_factory=list)
    lock: threading.RLock = field(default_factory=threading.RLock, repr=False)

    def __post_init__(self) -> None:
        self.workspace_root = self.workspace_root.expanduser().resolve()
        self.launch_root = self.launch_root.expanduser().resolve()
        self.workspace_root.mkdir(parents=True, exist_ok=True)
        self.launch_root.mkdir(parents=True, exist_ok=True)
        if self.active_root is not None:
            self.active_root = self.active_root.expanduser().resolve()

    def data_root(self) -> Path:
        with self.lock:
            return self.active_root or self.workspace_root

    def _refresh_process(self) -> None:
        if self.process is None:
            return
        returncode = self.process.poll()
        if returncode is None:
            return
        self.process_returncode = int(returncode)
        if self.process_finished_at is None:
            self.process_finished_at = _utc_now()
        if self.process_returncode == 0 and self.completion_error is None:
            self.completion_error = self._completed_run_error()
            if self.completion_error:
                self._record_completion_failure(self.completion_error)

    def _completed_run_error(self) -> str | None:
        """Reject a false-success child that produced no simulation results."""
        root = self.active_root
        if root is None or not self.command or "explore" not in self.command:
            return None

        setup_dir = root / "setup"
        required_setup = ("system.xml", "state.xml", "topology.pdb")
        missing_setup = [name for name in required_setup if not (setup_dir / name).is_file()]
        if missing_setup:
            setup_manifest = _json_mapping(setup_dir / "setup_parameters.json")
            detail = _manifest_note(setup_manifest)
            suffix = f" Details: {detail}" if detail else ""
            return (
                "Setup did not produce the files required for simulation "
                f"({', '.join(missing_setup)}).{suffix} See "
                f"{self.log_path or root / 'dashboard_launcher.log'} for the full log."
            )

        simulation_dir = root / "simulation"
        simulation_manifest = _json_mapping(
            simulation_dir / "simulation_parameters.json"
        )
        final_state = simulation_dir / "state_final.xml"
        simulated = (
            _is_nonempty_file(final_state)
            and simulation_manifest.get("platform_used") not in (None, "")
            and simulation_manifest.get("duration_ns_actual") is not None
        )
        if not simulated:
            detail = _manifest_note(simulation_manifest)
            suffix = f" Details: {detail}" if detail else ""
            return (
                "The workflow exited without producing a completed molecular "
                f"dynamics simulation.{suffix} See "
                f"{self.log_path or root / 'dashboard_launcher.log'} for the full log."
            )
        return None

    def _record_completion_failure(self, detail: str) -> None:
        """Make the overview and health panels reflect post-run validation."""
        if self.active_root is None:
            return
        try:
            from fastmdxplora.live.telemetry import TelemetryWriter, read_status

            status = read_status(self.active_root)
            states = status.get("stage_states")
            states = dict(states) if isinstance(states, dict) else {}
            setup_ready = all(
                (self.active_root / "setup" / name).is_file()
                for name in ("system.xml", "state.xml", "topology.pdb")
            )
            if setup_ready:
                stage = "production"
                states["production"] = "failed"
                states["analysis"] = "skipped"
                states["report"] = "skipped"
            else:
                stage = "setup"
                states["setup"] = "failed"
                for name in (
                    "minimization",
                    "nvt",
                    "npt",
                    "production",
                    "analysis",
                    "report",
                ):
                    states[name] = "skipped"
            writer = TelemetryWriter(self.active_root / "simulation", enabled=True)
            writer.write_status(
                stage=stage,
                status="failed",
                latest_error=detail,
                stage_states=states,
            )
            writer.event(detail, level="error")
        except Exception:  # noqa: BLE001 - status reporting must not mask the result
            return

    def snapshot(self) -> dict[str, Any]:
        with self.lock:
            self._refresh_process()
            running = self.process is not None and self.process.poll() is None
            status = "idle"
            if running:
                status = "running"
            elif self.completion_error:
                status = "failed"
            elif self.process is not None and self.process_returncode == 0:
                status = "completed"
            elif self.process is not None and self.process_returncode is not None:
                status = "failed"
            return {
                "mode": "home" if self.active_root is None else "run",
                "status": status,
                "active_run": str(self.active_root) if self.active_root else None,
                "workspace": str(self.workspace_root),
                "launch_root": str(self.launch_root),
                "process_running": running,
                "returncode": self.process_returncode,
                "error": self.completion_error,
                "started_at": self.process_started_at,
                "finished_at": self.process_finished_at,
                "log_path": str(self.log_path) if self.log_path else None,
                "command": list(self.command),
                "can_launch": not running,
            }

    def launch(self, payload: Mapping[str, Any], *, dashboard_url: str | None = None) -> dict[str, Any]:
        result = validate_launcher_payload(payload)
        if not result["valid"]:
            return result
        config = result["config"]
        with self.lock:
            self._refresh_process()
            if self.process is not None and self.process.poll() is None:
                return {
                    **result,
                    "valid": False,
                    "errors": {"run": "A FastMDXplora workflow is already running."},
                }

            environment_error = launcher_environment_error(config)
            if environment_error:
                return {
                    **result,
                    "valid": False,
                    "error": environment_error,
                    "errors": {"run": environment_error},
                }

            run_name = _slug(config["run_name"])
            output_dir = (self.launch_root / run_name).resolve()
            try:
                output_dir.relative_to(self.launch_root)
            except ValueError:
                return {
                    **result,
                    "valid": False,
                    "errors": {"run_name": "Run name resolves outside the launch directory."},
                }
            if output_dir.exists() and any(output_dir.iterdir()):
                return {
                    **result,
                    "valid": False,
                    "errors": {"run_name": f"Output folder already exists and is not empty: {output_dir}"},
                }
            output_dir.mkdir(parents=True, exist_ok=True)
            log_path = output_dir / "dashboard_launcher.log"
            command = build_launcher_command(config, output_dir)
            env = os.environ.copy()
            env["FASTMDX_DASHBOARD_ACTIVE"] = "1"
            env["FASTMDX_DASHBOARD_OUTPUT"] = str(output_dir)
            if dashboard_url:
                env["FASTMDX_DASHBOARD_URL"] = dashboard_url
            log_handle = log_path.open("a", encoding="utf-8", buffering=1)
            try:
                process = subprocess.Popen(
                    command,
                    cwd=str(self.launch_root),
                    env=env,
                    stdout=log_handle,
                    stderr=subprocess.STDOUT,
                    stdin=subprocess.DEVNULL,
                    shell=False,
                )
            except Exception:
                log_handle.close()
                raise
            # Popen owns an inherited OS handle. Closing our copy avoids a
            # long-lived Python file object while the child continues writing.
            log_handle.close()
            self.active_root = output_dir
            self.process = process
            self.process_started_at = _utc_now()
            self.process_finished_at = None
            self.process_returncode = None
            self.completion_error = None
            self.log_path = log_path
            self.command = command
            return {
                **result,
                "launched": True,
                "output": str(output_dir),
                "pid": process.pid,
                "command": command,
                "state": self.snapshot(),
            }

    def stop(self) -> dict[str, Any]:
        with self.lock:
            self._refresh_process()
            if self.process is None or self.process.poll() is not None:
                return {"stopped": False, "detail": "No workflow is currently running.", "state": self.snapshot()}
            self.process.terminate()
            return {"stopped": True, "detail": "Termination requested.", "state": self.snapshot()}
