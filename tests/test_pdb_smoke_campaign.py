"""Tests for the PDB smoke campaign runner."""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "run_pdb_smoke_campaign.py"
SPEC = importlib.util.spec_from_file_location("run_pdb_smoke_campaign", SCRIPT)
campaign = importlib.util.module_from_spec(SPEC)
assert SPEC and SPEC.loader
sys.modules[SPEC.name] = campaign
SPEC.loader.exec_module(campaign)


def _args(tmp_path: Path, **overrides):
    values = {
        "output_root": tmp_path / "campaign",
        "preset": "gentle",
        "nvt_steps": 1000,
        "npt_steps": 0,
        "production_steps": 1000,
        "trajectory_interval_steps": 100,
        "platform": "auto",
        "max_input_mb": 10.0,
        "max_setup_atoms": 0,
        "no_report": False,
        "verbose": False,
    }
    values.update(overrides)
    return SimpleNamespace(**values)


def _phase(name: str, status: str = "ok", message: str = ""):
    return SimpleNamespace(name=name, status=status, message=message)


class FakeFastMDXploraSuccess:
    """FastMDXplora stand-in that writes the artifacts validation expects."""

    captured_system = None
    captured_options = None

    def __init__(self, *, system, output_dir, verbose=False):
        self.system = system
        self.output_dir = Path(output_dir)
        self.verbose = verbose
        self.output_dir.mkdir(parents=True, exist_ok=True)
        FakeFastMDXploraSuccess.captured_system = system

    def explore(self, *, options, report=True):
        FakeFastMDXploraSuccess.captured_options = options
        setup = self.output_dir / "setup"
        sim = self.output_dir / "simulation"
        analysis = self.output_dir / "analysis"
        report_dir = self.output_dir / "report"
        for path in (setup, sim, analysis, report_dir):
            path.mkdir()
        (setup / "system.xml").write_text("<System><Particle/></System>", encoding="utf-8")
        (setup / "state.xml").write_text("<State/>", encoding="utf-8")
        (setup / "topology.pdb").write_text(
            # A tripeptide rather than a lone residue: one amino acid is
        # simultaneously N- and C-terminal, which AMBER has no template for.
        # These fixtures test pipeline mechanics, not force-field edge cases.
                "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N\n"
                "ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00  0.00           C\n"
        "ATOM      3  C   ALA A   1       2.009   1.420   0.000  1.00  0.00           C\n"
        "ATOM      4  O   ALA A   1       1.251   2.390   0.000  1.00  0.00           O\n"
        "ATOM      5  CB  ALA A   1       1.988  -0.773  -1.199  1.00  0.00           C\n"
        "ATOM      6  N   GLY A   2       3.332   1.549   0.000  1.00  0.00           N\n"
        "ATOM      7  CA  GLY A   2       3.972   2.849   0.000  1.00  0.00           C\n"
        "ATOM      8  C   GLY A   2       5.486   2.705   0.000  1.00  0.00           C\n"
        "ATOM      9  O   GLY A   2       6.008   1.593   0.000  1.00  0.00           O\n"
        "ATOM     10  N   ALA A   3       6.171   3.845   0.000  1.00  0.00           N\n"
        "ATOM     11  CA  ALA A   3       7.623   3.845   0.000  1.00  0.00           C\n"
        "ATOM     12  C   ALA A   3       8.174   5.265   0.000  1.00  0.00           C\n"
        "ATOM     13  O   ALA A   3       7.416   6.235   0.000  1.00  0.00           O\n"
        "ATOM     14  CB  ALA A   3       8.153   3.072  -1.199  1.00  0.00           C\n"
        "ATOM     15  OXT ALA A   3       9.400   5.400   0.000  1.00  0.00           O\n"
        "END\n",
            encoding="utf-8",
        )
        (setup / "setup_parameters.json").write_text(
            json.dumps({"notes": []}), encoding="utf-8"
        )
        (sim / "state_minimized.xml").write_text("<State/>", encoding="utf-8")
        (sim / "production.dcd").write_text("stub", encoding="utf-8")
        (sim / "simulation_parameters.json").write_text(
            json.dumps({"notes": []}), encoding="utf-8"
        )
        (analysis / "analysis_manifest.json").write_text("{}", encoding="utf-8")
        (report_dir / "report.md").write_text("# Report\n", encoding="utf-8")
        phases = [_phase("setup"), _phase("simulation"), _phase("analysis"), _phase("report")]
        manifest = {"phases": [{"name": p.name, "status": p.status} for p in phases]}
        (self.output_dir / "manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
        return [SimpleNamespace(output_dir=self.output_dir, phases=phases)]


class FakeFastMDXploraFailure(FakeFastMDXploraSuccess):
    def explore(self, *, options, report=True):
        setup = self.output_dir / "setup"
        setup.mkdir()
        (setup / "setup_parameters.json").write_text(
            json.dumps({"notes": ["PDBFixer unavailable: install via conda-forge"]}),
            encoding="utf-8",
        )
        phases = [_phase("setup", "ok"), _phase("simulation", "error", "OpenMM NaN")]
        manifest = {"phases": [{"name": p.name, "status": p.status} for p in phases]}
        (self.output_dir / "manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
        return [SimpleNamespace(output_dir=self.output_dir, phases=phases)]


def test_error_classification_handles_expected_categories() -> None:
    assert campaign.classify_failure("Temporary failure in name resolution") == (
        "DNS/download failure"
    )
    assert campaign.classify_failure("No template found for residue MSE") == (
        "unsupported residue/template failure"
    )
    assert campaign.classify_failure("Particle coordinate is NaN") == "OpenMM NaN"
    assert campaign.classify_failure("slides.pptx could not be written") == (
        "report generation failure"
    )
    assert campaign.classify_failure("Could not classify system input 'missing.pdb'") == (
        "missing input file"
    )


def test_summary_csv_and_json_are_written(tmp_path: Path) -> None:
    rows = [
        {
            "pdb_id": "1L2Y",
            "input_path": "1L2Y",
            "status": "ok",
            "output_dir": str(tmp_path / "campaign" / "1L2Y"),
        }
    ]
    campaign.write_summaries(rows, tmp_path / "campaign")
    assert (tmp_path / "campaign" / "campaign_summary.csv").exists()
    data = json.loads((tmp_path / "campaign" / "campaign_summary.json").read_text())
    assert data[0]["pdb_id"] == "1L2Y"


def test_run_one_uses_local_pdb_path_and_gentle_by_default(tmp_path: Path) -> None:
    local_pdb = tmp_path / "local.pdb"
    local_pdb.write_text("END\n", encoding="utf-8")
    with patch("fastmdxplora.FastMDXplora", FakeFastMDXploraSuccess), patch.object(
        campaign, "check_trajectory_coordinates", return_value=[]
    ):
        row = campaign.run_one(str(local_pdb), _args(tmp_path))

    assert row["status"] == "ok"
    assert row["input_path"] == str(local_pdb.resolve())
    assert FakeFastMDXploraSuccess.captured_system == str(local_pdb)
    assert FakeFastMDXploraSuccess.captured_options["simulation"]["preset"] == "gentle"


def test_missing_file_is_classified_without_crashing(tmp_path: Path) -> None:
    row = campaign._base_row(str(tmp_path / "missing.pdb"), tmp_path / "out")
    row.update(
        {
            "status": "failed",
            "failure_category": campaign.classify_failure("No such file or directory"),
            "error_message": "No such file or directory",
        }
    )
    assert row["failure_category"] == "missing input file"
    assert campaign.classify_bug_likelihood(row) == "expected limitation/input issue"


def test_campaign_runner_continues_and_summarizes_failure(tmp_path: Path) -> None:
    with patch("fastmdxplora.FastMDXplora", FakeFastMDXploraFailure), patch.object(
        campaign, "check_trajectory_coordinates", return_value=[]
    ):
        row = campaign.run_one("1L2Y", _args(tmp_path))

    assert row["status"] == "failed"
    assert row["failed_phase"] == "simulation"
    assert row["failure_category"] == "OpenMM NaN"
    assert row["bug_classification"] == "likely code bug"


def test_dependency_note_takes_precedence_over_missing_artifacts(tmp_path: Path) -> None:
    out = tmp_path / "run"
    setup = out / "setup"
    setup.mkdir(parents=True)
    (out / "manifest.json").write_text(
        json.dumps({"phases": [{"name": "setup", "status": "ok"}]}),
        encoding="utf-8",
    )
    (setup / "setup_parameters.json").write_text(
        json.dumps({"notes": ["PDBFixer unavailable: install via conda-forge"]}),
        encoding="utf-8",
    )
    validation = campaign.validate_output_dir(out)
    row = campaign._summarize_run_result(
        SimpleNamespace(output_dir=out, phases=[_phase("setup")]),
        validation,
    )
    assert row["status"] == "expected_limitation"
    assert row["failure_category"] == "missing dependency"
    assert row["bug_classification"] == "expected limitation/input issue"


def test_large_local_input_is_skipped_before_running_pipeline(tmp_path: Path) -> None:
    local_pdb = tmp_path / "large.pdb"
    local_pdb.write_bytes(b"x" * 2048)
    row = campaign.run_one(str(local_pdb), _args(tmp_path, max_input_mb=0.001))
    assert row["status"] == "skipped"
    assert row["failure_category"] == "too large for smoke settings"
    assert row["bug_classification"] == "expected limitation/input issue"


def test_parse_input_list_prefers_clean_local_lines(tmp_path: Path) -> None:
    pdb = tmp_path / "local.pdb"
    pdb.write_text("END\n", encoding="utf-8")
    input_list = tmp_path / "pdb_list.txt"
    input_list.write_text(f"# comment\n1L2Y\n{pdb}  # local\n1L2Y\n", encoding="utf-8")
    args = SimpleNamespace(input_list=[str(input_list)], inputs=[])
    assert campaign.parse_inputs(args) == ["1L2Y", str(pdb)]


@pytest.mark.skipif(
    not bool(__import__("os").environ.get("FASTMDX_RUN_OPENMM_TESTS")),
    reason="set FASTMDX_RUN_OPENMM_TESTS=1 to run real OpenMM/PDBFixer smoke tests",
)
def test_real_tiny_pdb_ids_when_openmm_enabled(tmp_path: Path) -> None:
    pytest.importorskip("openmm")
    pytest.importorskip("pdbfixer")
    args = _args(tmp_path, output_root=tmp_path / "real_campaign", no_report=True)
    rows = [campaign.run_one(pdb_id, args) for pdb_id in ("1L2Y", "1CRN")]
    assert all(row["status"] in {"ok", "expected_limitation"} for row in rows)
    campaign.write_summaries(rows, Path(args.output_root))
    assert (Path(args.output_root) / "campaign_summary.csv").exists()
