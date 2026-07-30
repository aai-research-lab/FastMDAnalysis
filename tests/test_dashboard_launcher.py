from __future__ import annotations

import json
import urllib.error
import urllib.request
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from fastmdxplora.live.launcher import (
    DashboardRuntime,
    build_launcher_command,
    launcher_environment_error,
    launcher_defaults,
    validate_launcher_payload,
)
from fastmdxplora.live.server import start_test_server


def _payload() -> dict:
    return {
        "system": "1L2Y",
        "run_name": "trpcage_test",
        "setup": {
            "ph": 7.4,
            "forcefield": "charmm36",
            "water_model": "auto",
            "ion_concentration_M": 0.15,
            "solvent_padding_nm": 1.0,
        },
        "simulation": {
            "minimize": True,
            "nvt_steps": 1000,
            "npt_steps": 1000,
            "production_steps": 10000,
            "timestep_fs": 2.0,
            "temperature_K": 300,
            "friction_per_ps": 1.0,
            "integrator": "langevin_middle",
            "platform": "CPU",
            "precision": "mixed",
            "trajectory_interval_steps": 100,
            "checkpoint_interval_steps": 1000,
            "telemetry_interval": 100,
        },
        "workflow": {
            "run_analysis": True,
            "run_report": True,
            "analyses": ["rmsd", "rg"],
            "report_document": True,
            "report_slides": True,
            "report_bundle": True,
        },
    }


def test_launcher_defaults_are_backend_derived() -> None:
    defaults = launcher_defaults()
    assert defaults["setup"]["forcefield"] == "charmm36"
    assert defaults["simulation"]["nvt_steps"] == 250_000
    assert "CPU" in defaults["choices"]["platforms"]


def test_launcher_validation_computes_durations() -> None:
    result = validate_launcher_payload(_payload())
    assert result["valid"] is True
    assert result["summary"]["production_ns"] == 0.02
    assert result["summary"]["trajectory_frames"] == 100


def test_launcher_validation_rejects_bad_values() -> None:
    payload = _payload()
    payload["system"] = ""
    payload["simulation"]["production_steps"] = 0
    result = validate_launcher_payload(payload)
    assert result["valid"] is False
    assert "system" in result["errors"]
    assert "simulation.production_steps" in result["errors"]


def test_launcher_command_uses_module_entrypoint(tmp_path: Path) -> None:
    result = validate_launcher_payload(_payload())
    command = build_launcher_command(result["config"], tmp_path / "out")
    assert command[1:4] == ["-m", "fastmdxplora.cli.main", "explore"]
    assert "--simulate-live-telemetry" in command
    assert "--dashboard" not in command
    assert "--analyze-analyses" in command


def test_launcher_environment_preflight_is_workflow_aware() -> None:
    imported: list[str] = []

    def fake_import(name: str):
        imported.append(name)
        if name == "pdbfixer":
            raise ImportError("missing")
        return SimpleNamespace()

    with patch("fastmdxplora.live.launcher.importlib.import_module", side_effect=fake_import):
        detail = launcher_environment_error(_payload())

    assert detail is not None
    assert "PDBFixer" in detail
    assert "MDTraj" not in detail
    assert imported == ["openmm", "openmm.app", "pdbfixer", "mdtraj"]

    payload = _payload()
    payload["workflow"]["run_analysis"] = False
    imported.clear()
    with patch("fastmdxplora.live.launcher.importlib.import_module", side_effect=fake_import):
        launcher_environment_error(payload)
    assert "mdtraj" not in imported


def test_runtime_launches_without_shell(tmp_path: Path) -> None:
    runtime = DashboardRuntime(
        workspace_root=tmp_path / "workspace",
        launch_root=tmp_path / "runs",
    )
    fake_process = SimpleNamespace(pid=42, poll=lambda: None, terminate=lambda: None)
    with (
        patch("fastmdxplora.live.launcher.launcher_environment_error", return_value=None),
        patch("fastmdxplora.live.launcher.subprocess.Popen", return_value=fake_process) as popen,
    ):
        result = runtime.launch(_payload(), dashboard_url="http://127.0.0.1:8765")
    assert result["launched"] is True
    kwargs = popen.call_args.kwargs
    assert kwargs["shell"] is False
    assert kwargs["env"]["FASTMDX_DASHBOARD_ACTIVE"] == "1"
    assert runtime.data_root().name == "trpcage_test"


def test_runtime_refuses_launch_when_simulation_dependencies_are_missing(
    tmp_path: Path,
) -> None:
    runtime = DashboardRuntime(
        workspace_root=tmp_path / "workspace",
        launch_root=tmp_path / "runs",
    )
    detail = "Simulation dependencies are unavailable: OpenMM, PDBFixer."
    with (
        patch("fastmdxplora.live.launcher.launcher_environment_error", return_value=detail),
        patch("fastmdxplora.live.launcher.subprocess.Popen") as popen,
    ):
        result = runtime.launch(_payload())
    assert result["valid"] is False
    assert result["error"] == detail
    assert result["errors"]["run"] == detail
    popen.assert_not_called()
    assert not (runtime.launch_root / "trpcage_test").exists()


def test_runtime_rejects_zero_exit_without_simulation_outputs(tmp_path: Path) -> None:
    runtime = DashboardRuntime(
        workspace_root=tmp_path / "workspace",
        launch_root=tmp_path / "runs",
    )
    run_root = runtime.launch_root / "incomplete"
    (run_root / "setup").mkdir(parents=True)
    (run_root / "simulation").mkdir()
    (run_root / "setup" / "setup_parameters.json").write_text(
        json.dumps({"notes": ["PDBFixer unavailable"]}),
        encoding="utf-8",
    )
    runtime.active_root = run_root
    runtime.log_path = run_root / "dashboard_launcher.log"
    runtime.command = ["python", "-m", "fastmdxplora.cli.main", "explore"]
    runtime.process = SimpleNamespace(poll=lambda: 0)

    state = runtime.snapshot()

    assert state["status"] == "failed"
    assert state["returncode"] == 0
    assert "Setup did not produce" in state["error"]
    live_status = json.loads(
        (run_root / "simulation" / "live_status.json").read_text(encoding="utf-8")
    )
    assert live_status["status"] == "failed"
    assert live_status["stage_states"]["setup"] == "failed"
    assert live_status["stage_states"]["production"] == "skipped"
    assert live_status["stage_states"]["analysis"] == "skipped"
    assert live_status["stage_states"]["report"] == "skipped"


def test_runtime_accepts_zero_exit_with_completed_simulation(tmp_path: Path) -> None:
    runtime = DashboardRuntime(
        workspace_root=tmp_path / "workspace",
        launch_root=tmp_path / "runs",
    )
    run_root = runtime.launch_root / "complete"
    setup_dir = run_root / "setup"
    simulation_dir = run_root / "simulation"
    setup_dir.mkdir(parents=True)
    simulation_dir.mkdir()
    for name in ("system.xml", "state.xml", "topology.pdb"):
        (setup_dir / name).write_text("ready", encoding="utf-8")
    (simulation_dir / "state_final.xml").write_text("<State />", encoding="utf-8")
    (simulation_dir / "simulation_parameters.json").write_text(
        json.dumps({"platform_used": "CPU", "duration_ns_actual": 0.02}),
        encoding="utf-8",
    )
    runtime.active_root = run_root
    runtime.command = ["python", "-m", "fastmdxplora.cli.main", "explore"]
    runtime.process = SimpleNamespace(poll=lambda: 0)

    state = runtime.snapshot()

    assert state["status"] == "completed"
    assert state["error"] is None


def test_home_server_exposes_launcher_apis(tmp_path: Path) -> None:
    server, url = start_test_server(tmp_path / "workspace", home_mode=True)
    try:
        with urllib.request.urlopen(url + "/") as response:
            html = response.read().decode("utf-8")
        assert "New Simulation" in html
        assert "/static/simulation-builder.js" in html
        with urllib.request.urlopen(url + "/api/app-state") as response:
            state = json.load(response)
        assert state["active_run"] is None
        with urllib.request.urlopen(url + "/api/launcher/defaults") as response:
            defaults = json.load(response)
        assert defaults["simulation"]["temperature_K"] == 300.0

        request = urllib.request.Request(
            url + "/api/launcher/validate",
            data=json.dumps(_payload()).encode("utf-8"),
            headers={"Content-Type": "application/json"},
            method="POST",
        )
        with urllib.request.urlopen(request) as response:
            validated = json.load(response)
        assert validated["valid"] is True
        assert "command" in validated
    finally:
        server.shutdown()
        server.server_close()
