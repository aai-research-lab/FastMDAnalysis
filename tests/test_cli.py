"""CLI smoke tests."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import pytest

from fastmdxplora.orchestrator import FastMDXplora
from fastmdxplora.cli.main import main


def _make_pdb_stub(tmp_path: Path) -> Path:
    p = tmp_path / "stub.pdb"
    p.write_text("ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N\n")
    return p


def test_cli_version_subprocess() -> None:
    result = subprocess.run(
        [sys.executable, "-m", "fastmdxplora.cli.main", "--version"],
        capture_output=True,
        text=True,
        encoding="utf-8",
    )
    assert result.returncode == 0
    assert "fastmdx" in result.stdout
    assert "FastMDXplora" in result.stdout


def test_cli_help_subprocess() -> None:
    result = subprocess.run(
        [sys.executable, "-m", "fastmdxplora.cli.main", "--help"],
        capture_output=True,
        text=True,
        encoding="utf-8",
    )
    assert result.returncode == 0
    assert "FastMDXplora" in result.stdout
    assert "explore" in result.stdout
    assert "xplore" in result.stdout
    assert "setup" in result.stdout
    assert "simulate" in result.stdout
    assert "analyze" in result.stdout
    assert "report" in result.stdout


def test_cli_no_args_starts_dashboard_home() -> None:
    cli_main_module = sys.modules["fastmdxplora.cli.main"]
    with patch.object(cli_main_module, "_cmd_dashboard_home", return_value=0) as home:
        rc = main([])
    assert rc == 0
    home.assert_called_once_with()


def test_cli_cite() -> None:
    rc = main(["--cite"])
    assert rc == 0


def test_cli_info() -> None:
    rc = main(["info"])
    assert rc == 0


def test_repo_root_launcher_runs_without_installation() -> None:
    repo_root = Path(__file__).resolve().parents[1]
    result = subprocess.run(
        [sys.executable, str(repo_root / "fastmdx"), "--version"],
        capture_output=True,
        text=True,
        encoding="utf-8",
        cwd=repo_root,
    )
    assert result.returncode == 0
    assert "FastMDXplora" in result.stdout


def test_cli_health_subcommand_invokes_health() -> None:
    with patch("health.main", return_value=0) as health_main:
        rc = main(["health"])
    assert rc == 0
    health_main.assert_called_once_with([])


def test_cli_install_subcommand_uses_local_checkout() -> None:
    with patch("fastmdxplora.install.bootstrap_environment") as bootstrap_mock:
        rc = main(["install"])
    assert rc == 0
    bootstrap_mock.assert_called_once_with(
        env_name="fastmdxplora",
        python_version="3.10",
        yes=False,
        force=False,
        package_name=".",
        editable=False,
    )


def test_cli_bootstrap_alias_uses_local_checkout_when_running_from_repo() -> None:
    with patch("fastmdxplora.install.bootstrap_environment") as bootstrap_mock:
        rc = main(["bootstrap"])
    assert rc == 0
    bootstrap_mock.assert_called_once_with(
        env_name="fastmdxplora",
        python_version="3.10",
        yes=False,
        force=False,
        package_name=".",
        editable=False,
    )


def test_cli_install_e_subcommand_uses_editable_checkout() -> None:
    with patch("fastmdxplora.install.bootstrap_environment") as bootstrap_mock:
        rc = main(["install-e"])
    assert rc == 0
    bootstrap_mock.assert_called_once_with(
        env_name="fastmdxplora",
        python_version="3.10",
        yes=False,
        force=False,
        package_name=".",
        editable=True,
    )


def test_cli_explore(tmp_path: Path) -> None:
    pdb = _make_pdb_stub(tmp_path)
    out = tmp_path / "run"
    rc = main(["explore", "-system", str(pdb), "--output", str(out), "--simulate-nvt-steps", "2", "--simulate-npt-steps", "2", "--simulate-production-steps", "4", "--simulate-trajectory-interval-steps", "1"])
    assert rc == 0
    assert (out / "manifest.json").exists()
    assert (out / "setup" / "setup_parameters.json").exists()
    assert (out / "simulation" / "simulation_parameters.json").exists()
    assert (out / "analysis" / "analysis_manifest.json").exists()
    assert (out / "report" / "report.md").exists()


def test_cli_xplore_is_alias(tmp_path: Path) -> None:
    pdb = _make_pdb_stub(tmp_path)
    out = tmp_path / "run"
    rc = main(["xplore", "-system", str(pdb), "--output", str(out), "--simulate-nvt-steps", "2", "--simulate-npt-steps", "2", "--simulate-production-steps", "4", "--simulate-trajectory-interval-steps", "1"])
    assert rc == 0


def test_cli_explore_with_pdb_id(tmp_path: Path) -> None:
    out = tmp_path / "run"

    # This test verifies a 4-char PDB ID is fetched from RCSB and routed
    # through the full CLI pipeline (setup -> simulation -> analysis ->
    # report). We mock the MD engine: running a freshly-solvated real
    # protein through a 2-step equilibration is numerically unstable (NaN),
    # which is not what this test checks. The mock writes a real (tiny)
    # trajectory from the solvated topology so the downstream analysis and
    # report phases run on genuine data, exercising the real wiring.
    def _fake_run_simulation(*, topology_pdb, output_dir, **kwargs):
        import mdtraj as md
        from fastmdxplora.simulation.runner import SimulationResult

        traj = md.load(str(topology_pdb))
        # Stack a few identical frames so analyses that need >1 frame work.
        multi = md.join([traj] * 4)
        traj_path = output_dir / "production.dcd"
        multi.save_dcd(str(traj_path))

        # The real run_simulation writes a topology copy into the simulation
        # output dir; the analysis phase loads it from there.
        sim_topology = output_dir / "topology.pdb"
        traj[0].save_pdb(str(sim_topology))

        final_state = output_dir / "final_state.xml"
        final_state.write_text("<State/>\n", encoding="utf-8")
        energy_csv = output_dir / "energy.csv"
        energy_csv.write_text("step,potential\n0,0.0\n", encoding="utf-8")
        log_file = output_dir / "simulation.log"
        log_file.write_text("mock simulation\n", encoding="utf-8")

        return SimulationResult(
            trajectory=traj_path,
            topology=sim_topology,
            final_state=final_state,
            energy_csv=energy_csv,
            log_file=log_file,
            platform_used="CPU",
            n_production_frames=multi.n_frames,
            duration_ns_actual=0.0,
        )

    with patch(
        "fastmdxplora.simulation.runner.run_simulation",
        side_effect=_fake_run_simulation,
    ):
        rc = main(["explore", "--system", "1L2Y", "--output", str(out)])
    assert rc == 0
    assert (out / "setup" / "setup_parameters.json").exists()
    assert (out / "report" / "report.md").exists()


def test_cli_explore_requires_input(tmp_path: Path) -> None:
    # No -s/--system and no --config: explore reports the error and returns
    # exit code 2 (rather than raising), consistent with other CLI errors.
    rc = main(["explore", "--output", str(tmp_path / "run")])
    assert rc == 2


def test_cli_per_phase_setup(tmp_path: Path) -> None:
    pdb = _make_pdb_stub(tmp_path)
    out = tmp_path / "run"
    rc = main(["setup", "-system", str(pdb), "--output", str(out)])
    assert rc == 0
    assert (out / "setup" / "setup_parameters.json").exists()


def test_cli_report_can_rerun_from_existing_output_without_system(tmp_path: Path) -> None:
    pdb = _make_pdb_stub(tmp_path)
    out = tmp_path / "run"
    out.mkdir()
    (out / "manifest.json").write_text(
        json.dumps({"system": str(pdb)}),
        encoding="utf-8",
    )

    rc = main(["report", "--output", str(out), "--no-slides", "--no-bundle"])

    assert rc == 0
    assert (out / "report" / "report.md").exists()


def test_cli_explore_no_report(tmp_path: Path) -> None:
    pdb = _make_pdb_stub(tmp_path)
    out = tmp_path / "run"
    rc = main(["explore", "-system", str(pdb), "--output", str(out), "--no-report", "--simulate-nvt-steps", "2", "--simulate-npt-steps", "2", "--simulate-production-steps", "4", "--simulate-trajectory-interval-steps", "1"])
    assert rc == 0
    assert not (out / "report").exists() or not any((out / "report").iterdir())


def test_cli_explore_include_and_exclude_mutually_exclusive(tmp_path: Path) -> None:
    pdb = _make_pdb_stub(tmp_path)
    rc = main(
        [
            "explore",
            "-system",
            str(pdb),
            "--output",
            str(tmp_path / "run"),
            "--include",
            "setup",
            "--exclude",
            "report",
        ]
    )
    assert rc == 2


def test_cli_dashboard_flags_parse_on_workflow_commands() -> None:
    from fastmdxplora.cli.main import _build_parser

    parser = _build_parser()
    parser.parse_args(["explore", "--system", "1L2Y", "--dashboard"])
    parser.parse_args(["explore", "--system", "1L2Y", "--live-dashboard"])
    parser.parse_args(
        [
            "explore",
            "--system",
            "1L2Y",
            "--dashboard",
            "--dashboard-host",
            "127.0.0.1",
            "--dashboard-port",
            "8877",
            "--dashboard-stop-on-complete",
        ]
    )
    for command in ("setup", "simulate", "analyze", "report"):
        parser.parse_args([command, "--system", "1L2Y", "--dashboard"])


class _FakeDashboardSession:
    def __init__(
        self,
        *,
        url: str = "http://127.0.0.1:8765",
        requested_port: int = 8765,
        port: int = 8765,
    ) -> None:
        self.url = url
        self.requested_port = requested_port
        self.port = port
        self.port_was_changed = requested_port != port
        self.stopped = False
        self.waited = False

    def stop(self) -> None:
        self.stopped = True

    def wait_forever(self) -> None:
        self.waited = True


def test_cli_dashboard_starts_server_and_prints_url(tmp_path: Path, capsys) -> None:
    pdb = _make_pdb_stub(tmp_path)
    out = tmp_path / "run"
    out.mkdir()
    (out / "manifest.json").write_text(json.dumps({"system": str(pdb)}), encoding="utf-8")
    session = _FakeDashboardSession()

    with patch(
        "fastmdxplora.live.server.start_dashboard_session",
        return_value=session,
    ) as start:
        rc = main(
            [
                "report",
                "--output",
                str(out),
                "--no-slides",
                "--no-bundle",
                "--dashboard",
                "--dashboard-stop-on-complete",
            ]
        )

    assert rc == 0
    start.assert_called_once()
    assert start.call_args.kwargs["output"] == out
    assert start.call_args.kwargs["host"] == "127.0.0.1"
    assert start.call_args.kwargs["port"] == 8765
    assert session.stopped is True
    text = capsys.readouterr().out
    assert "Live dashboard running at: http://127.0.0.1:8765" in text
    assert f"Watching output folder: {out}" in text
    assert "Open this URL in your browser to monitor the run." in text


def test_cli_explore_dashboard_does_not_read_missing_output_dir(
    tmp_path: Path,
) -> None:
    pdb = _make_pdb_stub(tmp_path)
    out = tmp_path / "run"
    session = _FakeDashboardSession()

    def fake_explore(self, *, dry_run=False):
        self.output_dir = out
        return [SimpleNamespace(status="ok")]

    with patch(
        "fastmdxplora.live.server.start_dashboard_session",
        return_value=session,
    ) as start, patch.object(FastMDXplora, "explore", fake_explore):
        rc = main(
            [
                "explore",
                "--system",
                str(pdb),
                "--output",
                str(out),
                "--include",
                "setup",
                "simulation",
                "analysis",
                "report",
                "--simulate-preset",
                "gentle",
                "--dashboard",
                "--dashboard-stop-on-complete",
            ]
        )

    assert rc == 0
    assert start.call_args.kwargs["output"] == out.resolve()
    assert session.stopped is True


def test_cli_phase_dashboard_uses_cli_output_path(tmp_path: Path) -> None:
    pdb = _make_pdb_stub(tmp_path)
    out = tmp_path / "phase_run"
    session = _FakeDashboardSession()

    with patch(
        "fastmdxplora.live.server.start_dashboard_session",
        return_value=session,
    ) as start, patch.object(
        FastMDXplora,
        "setup",
        return_value=SimpleNamespace(status="ok"),
    ):
        rc = main(
            [
                "setup",
                "--system",
                str(pdb),
                "--output",
                str(out),
                "--dashboard",
                "--dashboard-stop-on-complete",
            ]
        )

    assert rc == 0
    assert start.call_args.kwargs["output"] == out.resolve()
    assert session.stopped is True


@pytest.mark.parametrize("command", ["simulate", "analyze", "report"])
def test_cli_dashboard_uses_cli_output_path_for_other_phases(
    tmp_path: Path,
    command: str,
) -> None:
    pdb = _make_pdb_stub(tmp_path)
    out = tmp_path / f"{command}_run"
    session = _FakeDashboardSession()

    method_name = {"simulate": "simulate", "analyze": "analyze", "report": "report"}[command]
    with patch(
        "fastmdxplora.live.server.start_dashboard_session",
        return_value=session,
    ) as start, patch.object(
        FastMDXplora,
        method_name,
        return_value=SimpleNamespace(status="ok"),
    ):
        rc = main(
            [
                command,
                "--system",
                str(pdb),
                "--output",
                str(out),
                "--dashboard",
                "--dashboard-stop-on-complete",
            ]
        )

    assert rc == 0
    assert start.call_args.kwargs["output"] == out.resolve()
    assert session.stopped is True


def test_cli_dashboard_prints_port_conflict_and_host_warning(
    tmp_path: Path,
    capsys,
) -> None:
    pdb = _make_pdb_stub(tmp_path)
    out = tmp_path / "run"
    out.mkdir()
    (out / "manifest.json").write_text(json.dumps({"system": str(pdb)}), encoding="utf-8")
    session = _FakeDashboardSession(
        url="http://0.0.0.0:8766",
        requested_port=8765,
        port=8766,
    )

    with patch(
        "fastmdxplora.live.server.start_dashboard_session",
        return_value=session,
    ):
        rc = main(
            [
                "report",
                "--output",
                str(out),
                "--no-slides",
                "--no-bundle",
                "--dashboard",
                "--dashboard-host",
                "0.0.0.0",
                "--dashboard-port",
                "8765",
                "--dashboard-stop-on-complete",
            ]
        )

    assert rc == 0
    text = capsys.readouterr().out
    assert "Live dashboard running at: http://0.0.0.0:8766" in text
    assert "Requested port 8765 was busy, so FastMDXplora used 8766." in text
    assert "Warning: dashboard is bound to 0.0.0.0" in text


def test_cli_dashboard_implies_live_telemetry_for_simulation(tmp_path: Path) -> None:
    pdb = _make_pdb_stub(tmp_path)
    out = tmp_path / "run"
    session = _FakeDashboardSession()

    with patch(
        "fastmdxplora.live.server.start_dashboard_session",
        return_value=session,
    ), patch.object(
        FastMDXplora,
        "simulate",
        return_value=SimpleNamespace(status="ok"),
    ) as simulate:
        rc = main(
            [
                "simulate",
                "--system",
                str(pdb),
                "--output",
                str(out),
                "--dashboard",
                "--dashboard-stop-on-complete",
            ]
        )

    assert rc == 0
    assert simulate.call_args.kwargs["live_telemetry"] is True
    assert session.stopped is True
