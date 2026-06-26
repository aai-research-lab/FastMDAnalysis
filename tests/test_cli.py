"""CLI smoke tests."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path
from unittest.mock import patch

import pytest

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


def test_cli_no_args_shows_help() -> None:
    rc = main([])
    assert rc == 0


def test_cli_cite() -> None:
    rc = main(["--cite"])
    assert rc == 0


def test_cli_info() -> None:
    rc = main(["info"])
    assert rc == 0


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
def test_cli_splash_path(capsys) -> None:
    rc = main(["splash", "--path", "--no-preview"])
    assert rc == 0
    out = capsys.readouterr().out
    assert "graphical_abstract.png" in out


def test_cli_splash_preview(capsys) -> None:
    rc = main(["splash"])
    assert rc == 0
    out = capsys.readouterr().out
    assert "FastMDXplora" in out
    assert "Graphical abstract" in out