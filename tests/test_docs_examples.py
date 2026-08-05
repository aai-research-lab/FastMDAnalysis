"""Lightweight checks for documented command examples.

These tests keep README/docs snippets from drifting away from the CLI surface.
They intentionally parse commands rather than executing MD-heavy workflows.
"""

from __future__ import annotations

import shlex
import subprocess
import sys

import pytest

from fastmdxplora.cli.main import _build_parser
from fastmdxplora.config import validate_config
from scripts.run_pdb_smoke_campaign import build_parser as build_campaign_parser


@pytest.mark.parametrize(
    "command",
    [
        "explore --system protein.pdb --dry-run",
        "xplore -s 1L2Y --dry-run",
        (
            "explore -s protein.pdb --setup-ph 7.4 "
            "--simulate-duration-ns 100 --simulate-platform CUDA --dry-run"
        ),
        (
            "explore -s protein.pdb --include setup simulation "
            "--simulate-preset gentle --simulate-platform CPU --dry-run"
        ),
        (
            "explore --system local_pdbs/1L2Y.pdb "
            "--output local_runs/trpcage_live_full "
            "--include setup simulation analysis report "
            "--simulate-preset gentle --dashboard"
        ),
        "setup --system protein.pdb --ph 6.5 --box-shape octahedron",
        (
            "simulate --system protein.pdb --output ./trpcage_study "
            "--duration-ns 50.0 --platform CUDA --dashboard"
        ),
        "analyze --output ./trpcage_study --analyses rmsd rg --selection 'name CA'",
        "report --output ./trpcage_study --no-slides --dashboard",
        "gui --output ./trpcage_study --port 8765",
        "init-config --minimal -o study.yml",
        "info",
    ],
)
def test_documented_fastmdx_commands_parse(command: str) -> None:
    parser = _build_parser()
    parser.parse_args(shlex.split(command))


def test_documented_campaign_command_parses() -> None:
    parser = build_campaign_parser()
    args = parser.parse_args(
        shlex.split(
            "--input-list examples/pdb_list.txt "
            "--output-root runs/pdb_smoke_starter "
            "--preset gentle --continue-on-error"
        )
    )

    assert args.preset == "gentle"
    assert args.continue_on_error is True


def test_documented_module_entrypoint_fallback_runs() -> None:
    result = subprocess.run(
        [sys.executable, "-m", "fastmdxplora.cli.main", "--version"],
        check=True,
        text=True,
        capture_output=True,
    )

    assert "fastmdx" in result.stdout


def test_documented_configuration_shape_validates() -> None:
    validate_config(
        {
            "systems": [{"id": "trpcage", "system": "1L2Y"}],
            "setup": {"ph": 7.0, "forcefield": "charmm36"},
            "simulation": {"preset": "gentle", "duration_ns": 10},
            "analysis": {"scope": "solute", "include": ["rmsd", "rmsf", "rg"]},
            "report": {
                "title": "My study",
                "slides": True,
                "comparison": True,
                "region_highlights": [
                    {
                        "label": "example flexible loop",
                        "start": 3,
                        "end": 7,
                        "color": "#4E79A7",
                    }
                ],
            },
        },
        require_systems=True,
    )


def test_the_metadynamics_examples_in_the_docs_actually_plan() -> None:
    """A documented example the software refuses is worse than none.

    `docs/simulations.md` was written before walls existed and showed a
    ligand_distance block with no bound -- which the software then began
    refusing, correctly, leaving the page telling people to run something that
    stops. The examples are parsed and planned here so that cannot happen
    quietly again.
    """
    import re
    from pathlib import Path

    import pytest

    pytest.importorskip("yaml")
    import mdtraj as md
    import yaml

    from fastmdxplora.simulation.metadynamics import plan_from_config

    page = Path(__file__).resolve().parents[1] / "docs" / "simulations.md"
    if not page.is_file():  # pragma: no cover - the page is optional
        return

    top = md.Topology()
    chain = top.add_chain()
    for index in range(240):
        residue = top.add_residue("ALA", chain, resSeq=index + 1)
        for name in ("N", "CA", "C", "O"):
            top.add_atom(name, md.element.carbon, residue)
    ligand_chain = top.add_chain()
    ligand = top.add_residue("BNZ", ligand_chain, resSeq=900)
    for index in range(6):
        top.add_atom(f"C{index}", md.element.carbon, ligand)

    blocks = re.findall(r"```yaml\n(simulation:\n  metadynamics:.*?)```",
                        page.read_text(encoding="utf-8"), re.S)
    assert blocks, "the page should show at least one worked example"

    for block in blocks:
        spec = yaml.safe_load(block)["simulation"]["metadynamics"]
        spec.setdefault("ligand_resname", "BNZ")
        # Raises if the example is one the software would refuse.
        plan_from_config(spec, top)
