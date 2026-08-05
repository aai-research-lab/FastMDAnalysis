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


def test_every_page_is_in_the_documentation_navigation() -> None:
    """A page nobody can reach from the contents is a page nobody reads.

    The navigation was one flat list of thirteen entries under a single
    caption, with the production guide ahead of the phases it describes and a
    design note missing from it entirely. It is grouped now, in the order a
    new user meets things.
    """
    from pathlib import Path

    docs = Path(__file__).resolve().parents[1] / "docs"
    if not docs.is_dir():  # pragma: no cover - docs are optional
        return

    index = (docs / "index.md").read_text(encoding="utf-8")
    pages = {p.stem for p in docs.glob("*.md")} - {"index"}
    missing = sorted(page for page in pages if page not in index)
    assert not missing, f"these pages are in no toctree: {missing}"


def test_the_navigation_is_grouped() -> None:
    """Thirteen entries under one caption is a list, not an order."""
    from pathlib import Path

    index = (Path(__file__).resolve().parents[1] / "docs" / "index.md")
    if not index.is_file():  # pragma: no cover
        return
    captions = [line for line in index.read_text(encoding="utf-8").splitlines()
                if line.startswith(":caption:")]
    assert len(captions) >= 4, (
        f"the documentation should be grouped; found {len(captions)} section(s)"
    )


def test_the_readme_points_at_every_important_page() -> None:
    """The README's documentation section is the map. A capability documented
    on a page nothing links to is one nobody finds."""
    from pathlib import Path

    repo = Path(__file__).resolve().parents[1]
    readme = (repo / "README.md").read_text(encoding="utf-8")
    for page in ("installation", "getting_started", "phases", "simulations",
                 "interactions", "cli_reference", "configuration", "gui"):
        assert f"{page}.html" in readme, f"the README does not link {page}"


def test_the_gui_is_called_the_gui() -> None:
    """It is the FastMDXplora GUI, not "a browser". A browser is what you open
    it in, and the two are worth telling apart in a document somebody reads to
    learn what the software is called."""
    from pathlib import Path

    repo = Path(__file__).resolve().parents[1]
    pages = [repo / "README.md"] + sorted((repo / "docs").glob("*.md"))

    wrong = []
    for page in pages:
        if not page.is_file():
            continue
        text = page.read_text(encoding="utf-8")
        for phrase in ("browser interface", "the browser dashboard",
                       "watch one in a browser"):
            if phrase in text:
                wrong.append(f"{page.name}: {phrase!r}")
    assert not wrong, f"these call the GUI something else: {wrong}"


def test_every_command_the_docs_name_exists() -> None:
    """`fastmdx health --no-fix` was documented in three pages and had never
    existed. A documented command that is not there sends somebody to a usage
    error on their first attempt at diagnosing something.
    """
    import re
    from pathlib import Path

    from fastmdxplora.cli.main import _build_parser

    parser = _build_parser()
    commands = set(parser._subparsers._group_actions[0].choices)

    repo = Path(__file__).resolve().parents[1]
    pages = [repo / "README.md"] + sorted((repo / "docs").glob("*.md"))

    missing = {}
    for page in pages:
        if not page.is_file():
            continue
        for match in re.finditer(r"^\s*(?:\$ )?fastmdx\s+([a-z][a-z-]+)",
                                 page.read_text(encoding="utf-8"), re.M):
            name = match.group(1)
            if name not in commands and not name.startswith("-"):
                missing.setdefault(page.name, set()).add(name)
    assert not missing, f"these commands are documented and do not exist: {missing}"


def test_every_fastmdx_flag_the_docs_name_exists() -> None:
    """`fastmdx gui --open-browser` was documented; the flag is
    `--no-browser`, and the default is the opposite of what the page implied.
    """
    import re
    from pathlib import Path

    from fastmdxplora.cli.main import _build_parser

    known = set()

    def collect(parser):
        for action in parser._actions:
            known.update(action.option_strings)
            if hasattr(action, "choices") and isinstance(action.choices, dict):
                for child in action.choices.values():
                    collect(child)

    collect(_build_parser())

    repo = Path(__file__).resolve().parents[1]
    pages = [repo / "README.md"] + sorted((repo / "docs").glob("*.md"))

    #: Flags belonging to other tools, which the docs legitimately mention.
    OTHERS = {
        "--upgrade", "--no-deps", "--format", "--query-gpu", "--delete",
        "--input-list", "--output-root", "--continue-on-error", "--no-cache-dir",
        "--force-reinstall", "--dry-run", "--editable", "--user", "--prefix",
        "--version", "--help", "--extra-index-url", "--index-url",
    }

    missing = {}
    for page in pages:
        if not page.is_file():
            continue
        text = page.read_text(encoding="utf-8")
        # Only flags on a line that mentions fastmdx, so pip's and rsync's
        # are not counted as ours.
        for line in text.splitlines():
            if "fastmdx" not in line and not line.strip().startswith("--"):
                continue
            for match in re.finditer(r"(?<![\w-])(--[a-z][a-z0-9-]{3,})(?![\w-])",
                                     line):
                flag = match.group(1)
                if flag not in known and flag not in OTHERS:
                    missing.setdefault(page.name, set()).add(flag)
    assert not missing, f"these flags are documented and do not exist: {missing}"


def test_the_phases_page_lists_every_analysis() -> None:
    """It is the page everything else links to for what gets measured, and it
    fell four features behind: omega dihedrals, MDS, the per-residue SASA
    average and the clustering seed were all added without it changing.

    A reference that is missing what was added is worse than one that is
    obviously old, because nothing about it looks wrong.
    """
    from pathlib import Path

    import fastmdxplora.analysis  # noqa: F401
    from fastmdxplora.analysis.describe import describe_all

    page = Path(__file__).resolve().parents[1] / "docs" / "phases.md"
    if not page.is_file():  # pragma: no cover
        return
    text = page.read_text(encoding="utf-8")

    missing = sorted(name for name in describe_all() if f"`{name}`" not in text)
    assert not missing, f"these analyses are not in the phases page: {missing}"


def test_an_analysis_describes_what_it_currently_does() -> None:
    """The summary comes from the module docstring and reaches the GUI, the
    report and the docs. Three of them still said what the analysis did before
    it was extended -- dihedrals claimed phi/psi after omega was added, dimred
    claimed three methods after MDS made four.
    """
    import fastmdxplora.analysis  # noqa: F401
    from fastmdxplora.analysis import dihedrals, dimred, sasa
    from fastmdxplora.analysis.describe import explain_analysis

    dihedral_text = (explain_analysis("dihedrals").get("detail") or "").lower()
    assert "omega" in dihedral_text, (
        "dihedrals computes omega and its description should say so"
    )

    dimred_text = (explain_analysis("dimred").get("detail") or "").lower()
    assert "mds" in dimred_text

    sasa_text = (explain_analysis("sasa").get("detail") or "").lower()
    assert "mean" in sasa_text or "average" in sasa_text

    # And the constants they describe.
    assert "omega" in dihedrals.VALID_ANGLES
    assert "mds" in dimred.VALID_METHODS
    assert "average_residue" in sasa.VALID_MODES
