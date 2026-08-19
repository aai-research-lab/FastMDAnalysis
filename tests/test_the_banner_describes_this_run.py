"""The banner announces only what the command actually does.

`fastmdx analyze` reads a finished trajectory. It was printing a SIMULATION
block anyway -- 1,000,000 production steps and a 500-step frame interval taken
from the schema defaults, belonging to no run. The gate existed but read only
--include/--exclude, and a phase subcommand says nothing through those, so the
"where nothing says, show everything" fallback applied.
"""
import io
import sys
from contextlib import redirect_stdout
from unittest import mock

import pytest

from fastmdxplora.utils.presenter import SessionPresenter


def banner_for(argv: list[str]) -> str:
    buf = io.StringIO()
    with mock.patch.object(sys, "argv", ["fastmdx"] + argv):
        with redirect_stdout(buf):
            SessionPresenter().banner(
                System="3PTB.pdb", Output="runs/x", Version="0.0.0")
    return buf.getvalue()


@pytest.mark.parametrize("subcommand", ["analyze", "report"])
def test_a_command_that_does_not_simulate_says_nothing_about_simulating(
        subcommand):
    text = banner_for([subcommand, "--output", "runs/x"])
    assert "SIMULATION" not in text
    assert "Production" not in text
    assert "DCD Frames" not in text


@pytest.mark.parametrize("subcommand", ["analyze", "report", "simulate"])
def test_a_command_that_does_not_build_a_system_says_nothing_about_setup(
        subcommand):
    assert "SETUP" not in banner_for([subcommand, "--output", "runs/x"])


def test_setup_alone_announces_setup_alone():
    text = banner_for(["setup", "-s", "3PTB.pdb"])
    assert "SETUP" in text
    assert "SIMULATION" not in text


def test_explore_announces_both():
    text = banner_for(["explore", "-c", "study.yml"])
    assert "SETUP" in text
    assert "SIMULATION" in text
