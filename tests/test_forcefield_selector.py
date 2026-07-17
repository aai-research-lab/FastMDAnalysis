"""Tests for the named force-field selector.

These verify the selector's *logic and plumbing* — name resolution, defaults,
validation errors, manifest recording, config schema, and CLI flags. They do
not run OpenMM: whether a resolved XML set actually templates a given residue
is real chemistry, exercised by the gated end-to-end setup tests. The XML
filenames in the registry were verified against OpenMM's bundled data.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from fastmdxplora.setup.forcefields import (
    DEFAULT_FORCEFIELD,
    ForceFieldUnavailableError,
    ForceFieldChoice,
    available_forcefields,
    require_forcefield_xmls,
    resolve_forcefield,
)
from fastmdxplora.setup.pipeline import run as setup_run

from unittest.mock import MagicMock, patch

_MINI_PDB = (
    "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N  \n"
    "ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00  0.00           C  \n"
    "ATOM      3  C   ALA A   1       2.009   1.420   0.000  1.00  0.00           C  \n"
    "ATOM      4  O   ALA A   1       1.251   2.390   0.000  1.00  0.00           O  \n"
    "ATOM      5  CB  ALA A   1       2.000  -0.700   1.200  1.00  0.00           C  \n"
    "END\n"
)


@pytest.fixture
def stub_orchestrator(tmp_path: Path):
    """Minimal FastMDXplora stand-in with a tiny PDB system on disk."""
    pdb = tmp_path / "mini.pdb"
    pdb.write_text(_MINI_PDB)
    orch = MagicMock()
    orch.system = str(pdb)
    orch.output_dir = tmp_path / "proj"
    orch.output_dir.mkdir()
    orch._presenter = None
    return orch


# ---------------------------------------------------------------------------
# Resolver
# ---------------------------------------------------------------------------
class TestResolveForceField:
    def test_default_is_charmm36(self):
        assert DEFAULT_FORCEFIELD == "charmm36"

    def test_none_resolves_to_default(self):
        assert resolve_forcefield(None).name == "charmm36"

    def test_charmm36_xmls(self):
        c = resolve_forcefield("charmm36")
        assert c.xmls == ("charmm36.xml", "charmm36/water.xml")
        assert c.supports_ligand is False

    def test_charmm36m_requires_the_2024_c36m_xmls(self):
        c = resolve_forcefield("charmm36m")
        assert c.xmls == ("charmm36_2024.xml", "charmm36_2024/water.xml")
        assert c.requires_exact_xml is True
        assert "charmm36.xml" not in c.xmls

    def test_charmm36m_missing_xml_fails_without_fallback(self, tmp_path, monkeypatch):
        from fastmdxplora.setup import forcefields

        monkeypatch.setattr(forcefields, "_openmm_data_dir", lambda: tmp_path)
        with pytest.raises(ForceFieldUnavailableError, match="will not substitute"):
            require_forcefield_xmls(resolve_forcefield("charmm36m"))

    def test_charmm36m_exact_xmls_pass_availability_check(self, tmp_path, monkeypatch):
        from fastmdxplora.setup import forcefields

        (tmp_path / "charmm36_2024.xml").write_text("<ForceField/>")
        (tmp_path / "charmm36_2024").mkdir()
        (tmp_path / "charmm36_2024" / "water.xml").write_text("<ForceField/>")
        monkeypatch.setattr(forcefields, "_openmm_data_dir", lambda: tmp_path)
        require_forcefield_xmls(resolve_forcefield("charmm36m"))

    def test_amber14_xmls_and_water(self):
        c = resolve_forcefield("amber14")
        assert c.xmls == ("amber14-all.xml", "amber14/tip3p.xml")
        assert c.water_model == "tip3p"

    def test_amber_fb15_xmls(self):
        c = resolve_forcefield("amber-fb15")
        assert c.xmls == ("amberfb15.xml", "tip3p.xml")

    def test_case_insensitive(self):
        assert resolve_forcefield("CHARMM36").name == "charmm36"
        assert resolve_forcefield("  Amber14  ").name == "amber14"

    def test_unknown_raises_with_choices(self):
        with pytest.raises(ValueError) as exc:
            resolve_forcefield("nonexistent-ff")
        msg = str(exc.value)
        # Lists the valid choices to guide the user.
        assert "charmm36" in msg and "amber14" in msg

    def test_available_is_sorted(self):
        names = available_forcefields()
        assert names == tuple(sorted(names))
        assert "charmm36" in names

    def test_returns_frozen_choice(self):
        c = resolve_forcefield("charmm36")
        assert isinstance(c, ForceFieldChoice)
        with pytest.raises(Exception):
            c.name = "other"  # frozen dataclass


# ---------------------------------------------------------------------------
# Pipeline validation + manifest
# ---------------------------------------------------------------------------
class TestForceFieldPlumbing:
    def test_named_and_raw_together_raises(self, stub_orchestrator):
        out_dir = stub_orchestrator.output_dir / "setup"
        out_dir.mkdir()
        with pytest.raises(ValueError, match="either .*forcefield.* or .*force_field"):
            setup_run(
                orchestrator=stub_orchestrator,
                output_dir=out_dir,
                forcefield="charmm36",
                force_field=["charmm36.xml", "charmm36/water.xml"],
            )

    def test_unknown_named_ff_raises(self, stub_orchestrator):
        out_dir = stub_orchestrator.output_dir / "setup"
        out_dir.mkdir()
        with pytest.raises(ValueError, match="Unknown force field"):
            setup_run(
                orchestrator=stub_orchestrator,
                output_dir=out_dir,
                forcefield="bogus-ff",
            )

    def test_charmm36m_unavailable_fails_before_scaffolded_setup(
        self, stub_orchestrator, tmp_path, monkeypatch
    ):
        from fastmdxplora.setup import forcefields

        monkeypatch.setattr(forcefields, "_openmm_data_dir", lambda: tmp_path)
        out_dir = stub_orchestrator.output_dir / "setup"
        out_dir.mkdir()
        with pytest.raises(ForceFieldUnavailableError, match="will not substitute"):
            setup_run(
                orchestrator=stub_orchestrator,
                output_dir=out_dir,
                forcefield="charmm36m",
            )

    def test_named_amber14_recorded_resolved(self, stub_orchestrator):
        out_dir = stub_orchestrator.output_dir / "setup"
        out_dir.mkdir()
        # Mock the chemistry: this test verifies the resolved force field is
        # *recorded* in the manifest, not that OpenMM can template a residue
        # (real parameterization is covered by the gated end-to-end tests).
        from fastmdxplora.setup import prepare as _prepare_mod
        with patch.object(_prepare_mod, "prepare_system", return_value={}):
            setup_run(
                orchestrator=stub_orchestrator,
                output_dir=out_dir,
                forcefield="amber14",
            )
        manifest = json.loads(
            (out_dir / "setup_parameters.json").read_text(encoding="utf-8"))
        resolved = manifest["resolved_forcefield"]
        assert resolved["source"] == "named"
        assert resolved["name"] == "amber14"
        assert resolved["xmls"] == ["amber14-all.xml", "amber14/tip3p.xml"]
        assert resolved["water_model"] == "tip3p"

    def test_raw_list_recorded_as_explicit(self, stub_orchestrator):
        out_dir = stub_orchestrator.output_dir / "setup"
        out_dir.mkdir()
        from fastmdxplora.setup import prepare as _prepare_mod
        with patch.object(_prepare_mod, "prepare_system", return_value={}):
            setup_run(
                orchestrator=stub_orchestrator,
                output_dir=out_dir,
                force_field=["amber14-all.xml", "amber14/tip3pfb.xml"],
            )
        manifest = json.loads(
            (out_dir / "setup_parameters.json").read_text(encoding="utf-8"))
        resolved = manifest["resolved_forcefield"]
        assert resolved["source"] == "explicit_xml_list"
        assert resolved["xmls"] == ["amber14-all.xml", "amber14/tip3pfb.xml"]


# ---------------------------------------------------------------------------
# Config schema + CLI
# ---------------------------------------------------------------------------
class TestForceFieldSchemaAndCLI:
    def test_schema_has_forcefield_field(self):
        from fastmdxplora.config.schema import PHASE_SCHEMAS
        setup = PHASE_SCHEMAS["setup"]
        f = setup.get("forcefield")
        assert f is not None
        assert f.default == "charmm36"

    def test_cli_has_setup_forcefield_flag(self):
        from fastmdxplora.cli.main import _SETUP_OPTIONS
        flags = {flag for flag, _kw, _opts in _SETUP_OPTIONS}
        assert "forcefield" in flags
        # choices are constrained to the registered names
        for flag, _kw, opts in _SETUP_OPTIONS:
            if flag == "forcefield":
                assert {"charmm36", "charmm36m", "amber14", "amber-fb15"} <= set(opts["choices"])
