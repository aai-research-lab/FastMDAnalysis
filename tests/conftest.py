"""Shared test isolation for optional chemistry backends."""

from __future__ import annotations

import importlib

import pytest


@pytest.fixture(autouse=True)
def allow_mocked_cli_pipelines_without_chemistry(monkeypatch: pytest.MonkeyPatch) -> None:
    """Let CLI wiring tests exercise mocked phases without OpenMM/PDBFixer.

    The CLI performs an early dependency preflight for real setup/simulation
    runs. These tests replace the simulation/setup work with mocks, so the
    preflight would incorrectly prevent them from testing downstream routing
    on CI images that intentionally omit the optional chemistry stack.
    """
    cli_main = importlib.import_module("fastmdxplora.cli.main")
    monkeypatch.setattr(cli_main, "_missing_chemistry_backends", lambda: [])
