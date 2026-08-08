"""Preparing each structure once, so the assertions can be many and cheap.

Setup takes seconds to minutes per structure and fetches from RCSB, so it
runs once per session per structure and every test reads the same output.
That also means a failure in preparation is reported once, against the
structure, rather than repeated against each thing that wanted to look at it.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path

import pytest

from .structures import CORPUS, Expectation


@dataclass
class Prepared:
    """What a setup run left behind, or why it did not."""

    expectation: Expectation
    directory: Path
    returncode: int
    output: str

    @property
    def succeeded(self) -> bool:
        return self.returncode == 0

    @property
    def prepared_pdb(self) -> Path:
        return self.directory / "setup" / "prepared.pdb"

    @property
    def solvated_pdb(self) -> Path:
        return self.directory / "setup" / "solvated.pdb"

    def refusal(self) -> str:
        """The message a failed run ended with, for tests about refusals."""
        return self.output


def _workspace() -> Path:
    """Where prepared structures go.

    Honours FASTMDXPLORA_VALIDATION_DIR so a developer can keep the outputs
    between runs and look at them, which is the point of the exercise: these
    tests exist because looking at the structure found what the suite did not.
    """
    override = os.environ.get("FASTMDXPLORA_VALIDATION_DIR")
    if override:
        path = Path(override).expanduser()
        path.mkdir(parents=True, exist_ok=True)
        return path
    return Path(tempfile.mkdtemp(prefix="fastmdx-validation-"))


@pytest.fixture(scope="session")
def workspace() -> Path:
    return _workspace()


def _run_setup(expectation: Expectation, workspace: Path) -> Prepared:
    directory = workspace / expectation.label.replace("/", "_")
    if directory.exists():
        shutil.rmtree(directory)

    command = [
        sys.executable, "-m", "fastmdxplora.cli.main", "explore",
        "--system", expectation.pdb_id,
        "--output", str(directory),
        "--include", "setup",
    ]
    if expectation.chains:
        command += ["--setup-chains", expectation.chains]
    if expectation.membrane:
        command += ["--setup-membrane", expectation.membrane,
                    "--setup-membrane-orient"]

    completed = subprocess.run(
        command, capture_output=True, text=True, timeout=1800)
    return Prepared(
        expectation=expectation,
        directory=directory,
        returncode=completed.returncode,
        output=completed.stdout + completed.stderr,
    )


@pytest.fixture(scope="session")
def prepared(workspace: Path):
    """Prepare a structure once, by kind, and remember the result."""
    cache: dict[str, Prepared] = {}

    def _for(kind: str) -> Prepared:
        if kind not in cache:
            expectation = next(e for e in CORPUS if e.kind == kind)
            cache[kind] = _run_setup(expectation, workspace)
        return cache[kind]

    return _for


def load(path: Path):
    """Read a structure, or skip the test saying which file was missing."""
    import mdtraj as md

    if not path.is_file():
        pytest.fail(f"{path.name} was not written")
    return md.load(str(path))
