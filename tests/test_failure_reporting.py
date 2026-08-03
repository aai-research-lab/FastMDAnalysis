"""A phase that could not do its job must say so through the exit code.

This file exists because the opposite was true four separate times: a
heterogen removed silently, a refusal whose reason was dropped by the
per-phase command, a refusal downgraded to a warning that still reported
success, and a failed download that produced an empty directory and exit 0.

Each was found by accident. The contract is stated here so the next one
fails a test instead.

The distinction the tests draw is between two things that look alike:

* an absent optional dependency, which is a choice the user made when
  installing, and degrades gracefully; and
* a failed operation, which is a failure and must be reported.
"""

from __future__ import annotations

import contextlib
import io
import logging
from pathlib import Path

import pytest

from fastmdxplora.cli.main import main

MINIMAL_PDB = (
    "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N\n"
    "ATOM      2  CA  ALA A   1       1.400   0.000   0.000  1.00  0.00           C\n"
    "END\n"
)


def _run(argv) -> tuple[int, str]:
    """Run a command, returning its exit code and everything it logged."""
    captured = io.StringIO()
    handler = logging.StreamHandler(captured)
    package_logger = logging.getLogger("fastmdx")
    package_logger.addHandler(handler)
    stdout = io.StringIO()
    try:
        with contextlib.redirect_stdout(stdout):
            code = main(argv)
    except SystemExit as exc:  # pragma: no cover - argparse paths
        code = int(exc.code or 0)
    finally:
        package_logger.removeHandler(handler)
    return code, captured.getvalue() + stdout.getvalue()


class TestFailuresAreReported:
    def test_an_unresolvable_identifier_fails(self, tmp_path) -> None:
        """A mistyped PDB code produced an empty directory and exit 0.

        The manifest recorded the failure, but nothing a caller reads did.
        """
        code, output = _run([
            "setup", "--system", "ZZZZ", "--output", str(tmp_path / "run"),
        ])
        assert code == 1
        assert "resolve input" in output.lower() or "failed" in output.lower()

    def test_a_missing_file_fails(self, tmp_path) -> None:
        code, _ = _run([
            "setup", "--system", str(tmp_path / "absent.pdb"),
            "--output", str(tmp_path / "run"),
        ])
        assert code == 1

    def test_the_manifest_survives_a_failure(self, tmp_path) -> None:
        """Reporting the failure must not discard what was learned first."""
        out = tmp_path / "run"
        _run(["setup", "--system", "ZZZZ", "--output", str(out)])
        assert (out / "setup" / "setup_parameters.json").is_file()


class TestGracefulDegradationIsPreserved:
    """Not installing an optional backend is a choice, not a failure."""

    def test_a_missing_chemistry_stack_does_not_fail_the_run(self, tmp_path) -> None:
        if _openmm_present():
            pytest.skip("only meaningful without the chemistry stack installed")

        import json

        structure = tmp_path / "in.pdb"
        structure.write_text(MINIMAL_PDB, encoding="utf-8")
        out = tmp_path / "run"
        code, _ = _run(["setup", "--system", str(structure), "--output", str(out)])

        assert code == 0, "an uninstalled optional backend is not a failed run"

        # The presenter is silent without a terminal, so the manifest is where
        # the reason is checked: the run must still say why it did less.
        parameters = json.loads(
            (out / "setup" / "setup_parameters.json").read_text(encoding="utf-8")
        )
        notes = " ".join(str(n) for n in parameters.get("notes", []))
        assert "openmm" in notes.lower() or "pdbfixer" in notes.lower(), (
            f"the manifest should record why setup stopped early, got: {notes!r}"
        )


def _openmm_present() -> bool:
    try:
        import openmm  # noqa: F401
    except ImportError:
        return False
    return True


class TestReasonsReachTheUser:
    """A refusal is only useful if its reason is visible."""

    def test_a_single_phase_reports_why_it_failed(self, tmp_path) -> None:
        """`explore` reported through the orchestrator; `setup` did not."""
        structure = tmp_path / "in.pdb"
        structure.write_text(
            MINIMAL_PDB.replace(
                "END\n",
                "HETATM   10  C1  BNZ A 200       5.000   5.000   5.000  1.00  0.00           C\nEND\n",
            ),
            encoding="utf-8",
        )
        code, output = _run([
            "setup", "--system", str(structure),
            "--output", str(tmp_path / "run"), "--heterogens", "auto",
        ])
        assert code == 1
        # The reason, not merely the fact.
        assert "PDB identifier" in output
