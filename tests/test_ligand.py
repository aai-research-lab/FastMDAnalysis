"""Tests for ligand / cofactor support (protein-ligand systems).

These verify the ligand feature's *logic and plumbing* — format detection,
the optional-dependency guard, force-field coherence validation, the
`amber-openff` registry entry, schema/CLI/manifest wiring. They do not run
OpenFF: real small-molecule parameterization (loading an SDF, building a
SystemGenerator system) is gated and verified end-to-end on a host with the
OpenFF stack installed. Plumbing tests mock the chemistry.
"""

from __future__ import annotations

import json
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from fastmdxplora.setup import ligand as ligand_mod
from fastmdxplora.setup import prepare as prepare_mod
from fastmdxplora.setup.forcefields import resolve_forcefield
from fastmdxplora.setup.ligand import (
    LigandError,
    detect_ligand_format,
)
from fastmdxplora.setup.pipeline import run as setup_run


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
    pdb = tmp_path / "mini.pdb"
    pdb.write_text(_MINI_PDB)
    orch = MagicMock()
    orch.system = str(pdb)
    orch.output_dir = tmp_path / "proj"
    orch.output_dir.mkdir()
    orch._presenter = None
    return orch


# ---------------------------------------------------------------------------
# Format detection
# ---------------------------------------------------------------------------
class TestDetectLigandFormat:
    def test_sdf(self):
        assert detect_ligand_format("lig.sdf") == "sdf"

    def test_mol2(self):
        assert detect_ligand_format("/path/to/lig.MOL2") == "mol2"

    def test_pdb_rejected(self):
        with pytest.raises(LigandError, match="Unsupported ligand format"):
            detect_ligand_format("lig.pdb")

    def test_no_extension_rejected(self):
        with pytest.raises(LigandError):
            detect_ligand_format("ligand")


# ---------------------------------------------------------------------------
# Net-charge inference (no OpenFF needed — uses a fake molecule)
# ---------------------------------------------------------------------------
class TestInferNetCharge:
    def _mol_with_charge(self, value):
        mol = MagicMock()
        mol.total_charge = MagicMock(magnitude=value)
        return mol

    def test_integer_charge(self):
        assert ligand_mod._infer_net_charge(self._mol_with_charge(-1.0)) == -1

    def test_zero_charge(self):
        assert ligand_mod._infer_net_charge(self._mol_with_charge(0.0)) == 0

    def test_non_integer_returns_none(self):
        assert ligand_mod._infer_net_charge(self._mol_with_charge(0.4)) is None

    def test_plain_number_charge(self):
        mol = MagicMock()
        mol.total_charge = 2  # no .magnitude attribute
        assert ligand_mod._infer_net_charge(mol) == 2


# ---------------------------------------------------------------------------
# Import guard
# ---------------------------------------------------------------------------
class TestOpenFFGuard:
    def test_missing_openff_raises_clear_error(self):
        # Simulate OpenFF not installed.
        import builtins
        real_import = builtins.__import__

        def fake_import(name, *a, **k):
            if name.startswith("openff"):
                raise ImportError("no openff")
            return real_import(name, *a, **k)

        with patch.object(builtins, "__import__", side_effect=fake_import):
            with pytest.raises(LigandError, match="OpenFF toolkit"):
                ligand_mod._import_openff()


# ---------------------------------------------------------------------------
# Force-field coherence (registry + validation)
# ---------------------------------------------------------------------------
class TestForceFieldCoherence:
    def test_amber_openff_registered_and_ligand_capable(self):
        c = resolve_forcefield("amber-openff")
        assert c.supports_ligand is True
        assert c.small_molecule_forcefield == "openff-2.2.1"
        # The full AMBER14 bundle rather than the protein file alone, so a
        # structure containing DNA or RNA can be prepared at all.
        assert c.xmls == ("amber14-all.xml", "amber14/tip3p.xml")

    def test_protein_only_ffs_not_ligand_capable(self):
        for name in ("charmm36", "amber14", "amber-fb15"):
            assert resolve_forcefield(name).supports_ligand is False
            assert resolve_forcefield(name).small_molecule_forcefield is None

    def test_ligand_with_non_ligand_ff_raises(self, stub_orchestrator, tmp_path):
        lig = tmp_path / "lig.sdf"
        lig.write_text("fake sdf\n")
        out_dir = stub_orchestrator.output_dir / "setup"
        out_dir.mkdir()
        # charmm36 doesn't support ligands -> clear error before any chemistry.
        with pytest.raises(ValueError, match="does not support ligands"):
            setup_run(
                orchestrator=stub_orchestrator,
                output_dir=out_dir,
                forcefield="charmm36",
                ligand=str(lig),
            )

    def test_ligand_with_raw_xml_list_raises(self, stub_orchestrator, tmp_path):
        lig = tmp_path / "lig.sdf"
        lig.write_text("fake sdf\n")
        out_dir = stub_orchestrator.output_dir / "setup"
        out_dir.mkdir()
        with pytest.raises(ValueError, match="named .*forcefield. selector"):
            setup_run(
                orchestrator=stub_orchestrator,
                output_dir=out_dir,
                force_field=["amber14/protein.ff14SB.xml", "amber14/tip3p.xml"],
                ligand=str(lig),
            )


# ---------------------------------------------------------------------------
# Manifest + schema + CLI plumbing (chemistry mocked)
# ---------------------------------------------------------------------------
class TestLigandPlumbing:
    def test_ligand_recorded_in_manifest(self, stub_orchestrator, tmp_path):
        lig = tmp_path / "drug.sdf"
        lig.write_text("fake sdf\n")
        out_dir = stub_orchestrator.output_dir / "setup"
        out_dir.mkdir()
        with patch.object(prepare_mod, "prepare_system", return_value={}):
            setup_run(
                orchestrator=stub_orchestrator,
                output_dir=out_dir,
                forcefield="amber-openff",
                ligand=str(lig),
                ligand_name="DRG",
                ligand_net_charge=-1,
            )
        manifest = json.loads(
            (out_dir / "setup_parameters.json").read_text(encoding="utf-8"))
        rec = manifest["resolved_forcefield"]
        assert rec["name"] == "amber-openff"
        assert rec["supports_ligand"] is True
        assert rec["ligand"]["name"] == "DRG"
        assert rec["ligand"]["net_charge"] == -1
        assert rec["ligand"]["small_molecule_forcefield"] == "openff-2.2.1"
        assert rec["ligand"]["files"] == [str(lig)]

    def test_schema_has_ligand_fields(self):
        from fastmdxplora.config.schema import PHASE_SCHEMAS
        setup = PHASE_SCHEMAS["setup"]
        for f in ("ligand", "ligand_forcefield", "ligand_name", "ligand_net_charge"):
            assert setup.get(f) is not None, f"missing schema field {f}"

    def test_cli_has_ligand_flags(self):
        from fastmdxplora.cli.main import _SETUP_OPTIONS, _build_parser

        flags = {flag for flag, _kw, _opts in _SETUP_OPTIONS}
        assert {"ligand", "ligand-forcefield", "ligand-name", "ligand-net-charge"} <= flags
        # The schema declares the accepted force fields; the parser reads them.
        setup = _build_parser()._subparsers._group_actions[0].choices["setup"]
        choices = next(a.choices for a in setup._actions if a.dest == "forcefield")
        assert "amber-openff" in choices

    def test_normalize_ligands(self):
        assert prepare_mod._normalize_ligands(None) == []
        assert prepare_mod._normalize_ligands("a.sdf") == ["a.sdf"]
        assert prepare_mod._normalize_ligands(["a.sdf", "b.sdf"]) == ["a.sdf", "b.sdf"]

    def test_add_ligand_sets_residue_name(self, monkeypatch):
        """_add_ligand_to_modeller must rename the ligand residue to the
        configured name. OpenFF's to_openmm() names it 'UNK'; downstream
        resname selection (scope, ligand_rmsd) needs the real name."""
        from unittest.mock import MagicMock

        # Fake openff.units.openmm.to_openmm
        import sys, types
        fake_units = types.ModuleType("openff.units.openmm")
        fake_units.to_openmm = lambda conf: conf
        monkeypatch.setitem(sys.modules, "openff.units.openmm", fake_units)
        # also ensure parent packages exist
        for mod in ("openff", "openff.units"):
            if mod not in sys.modules:
                monkeypatch.setitem(sys.modules, mod, types.ModuleType(mod))

        # Ligand topology with one 'UNK' residue
        lig_res = MagicMock(); lig_res.name = "UNK"
        omm_topology = MagicMock()
        omm_topology.residues.return_value = [lig_res]
        off_topology = MagicMock()
        off_topology.to_openmm.return_value = omm_topology
        ligand_mol = MagicMock()
        ligand_mol.to_topology.return_value = off_topology
        ligand_mol.conformers = [MagicMock()]

        # Modeller: 2 protein residues before add; after add, a 3rd (ligand).
        protein_residues = [MagicMock(name="r1"), MagicMock(name="r2")]
        merged_ligand_res = MagicMock(); merged_ligand_res.name = "UNK"
        modeller = MagicMock()
        modeller.topology.getNumResidues.return_value = 2
        modeller.topology.residues.return_value = protein_residues + [merged_ligand_res]

        prepare_mod._add_ligand_to_modeller(
            {}, modeller, ligand_mol, ligand_name="LIG"
        )
        # The input-topology residue was renamed...
        assert lig_res.name == "LIG"
        # ...and the merged-topology ligand residue was re-asserted to LIG.
        assert merged_ligand_res.name == "LIG"


# ---------------------------------------------------------------------------
# Ligand-protein clash detection (pure geometry)
# ---------------------------------------------------------------------------
class TestClashCheck:
    def _modeller_with(self, coords_nm):
        """Build a stub modeller whose positions return the given nm coords."""
        modeller = MagicMock()
        modeller.topology.getNumAtoms.return_value = len(coords_nm)
        # positions.value_in_unit(nm) -> list of (x,y,z) tuples
        quantity = MagicMock()
        quantity.value_in_unit.return_value = coords_nm
        modeller.positions = quantity
        return modeller

    def test_clash_raises(self):
        # 1 protein atom at origin, 1 ligand atom 0.05 nm away -> clash.
        coords = [(0.0, 0.0, 0.0), (0.05, 0.0, 0.0)]
        modeller = self._modeller_with(coords)
        with pytest.raises(ValueError, match="clashes with the protein"):
            prepare_mod._check_ligand_clashes(
                modeller, 1, MagicMock(),
                threshold_nm=0.15, ligand_name="LIG",
            )

    def test_no_clash_passes(self):
        # protein at origin, ligand 1.0 nm away -> fine.
        coords = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0)]
        modeller = self._modeller_with(coords)
        # Should not raise.
        prepare_mod._check_ligand_clashes(
            modeller, 1, MagicMock(),
            threshold_nm=0.15, ligand_name="LIG",
        )

    def test_threshold_boundary(self):
        # Exactly at threshold is allowed (strict less-than is a clash).
        coords = [(0.0, 0.0, 0.0), (0.15, 0.0, 0.0)]
        modeller = self._modeller_with(coords)
        prepare_mod._check_ligand_clashes(
            modeller, 1, MagicMock(),
            threshold_nm=0.15, ligand_name="LIG",
        )

    def test_mocked_positions_skip_gracefully(self):
        # If positions aren't real numeric quantities, the check is skipped.
        modeller = MagicMock()
        modeller.positions.value_in_unit.side_effect = TypeError("not real")
        # Should not raise.
        prepare_mod._check_ligand_clashes(
            modeller, 1, MagicMock(),
            threshold_nm=0.15, ligand_name="LIG",
        )

    def test_schema_and_cli_have_clash_options(self):
        from fastmdxplora.config.schema import PHASE_SCHEMAS
        from fastmdxplora.cli.main import _SETUP_OPTIONS
        setup = PHASE_SCHEMAS["setup"]
        assert setup.get("check_ligand_clashes") is not None
        assert setup.get("ligand_clash_threshold_nm") is not None
        flags = {flag for flag, _kw, _opts in _SETUP_OPTIONS}
        assert "no-ligand-clash-check" in flags
        assert "ligand-clash-threshold-nm" in flags


class TestHowTheLigandsChemistryWasKnown:
    """A trajectory carries coordinates. Whether a nitrogen donates, whether a
    ring is aromatic, whether a group is charged -- that is chemistry.

    Other profilers perceive it from the coordinates every time, because they
    accept arbitrary PDB files. This one often knows: setup resolves it from
    the Chemical Component Dictionary and settles the protonation against the
    pocket. Where it does not know, it says so, because a wrong bond order
    moves a hydrogen and invents or destroys a hydrogen bond.
    """

    @pytest.fixture(autouse=True)
    def _needs_rdkit(self):
        """A ligand's chemistry needs RDKit, which is the [ligand] extra.

        What is being tested is chemistry -- bond orders, aromaticity, formal
        charge -- so there is nothing useful to assert without it. A mock
        would only check that the mock was called.
        """
        pytest.importorskip("rdkit", reason="requires the [ligand] extra")


    @staticmethod
    def _ligand(smiles):
        import numpy as np
        import mdtraj as md
        from rdkit import Chem
        from rdkit.Chem import AllChem

        mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
        AllChem.EmbedMolecule(mol, randomSeed=3)
        elements = {"C": md.element.carbon, "O": md.element.oxygen,
                    "N": md.element.nitrogen, "H": md.element.hydrogen}
        top = md.Topology()
        chain = top.add_chain()
        residue = top.add_residue("LIG", chain, resSeq=1)
        for index, atom in enumerate(mol.GetAtoms()):
            top.add_atom(f"{atom.GetSymbol()}{index}",
                         elements[atom.GetSymbol()], residue)
        conf = mol.GetConformer()
        xyz = np.array(
            [[list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())]],
            dtype=np.float32) / 10.0
        traj = md.Trajectory(xyz=xyz, topology=top)
        return mol, traj, list(range(mol.GetNumAtoms()))

    @staticmethod
    def _double_bonds(mol):
        return sum(1 for b in mol.GetBonds() if b.GetBondTypeAsDouble() == 2)

    def test_a_supplied_file_is_used_and_named(self, tmp_path) -> None:
        from rdkit import Chem

        from fastmdxplora.analysis.ligand_chemistry import resolve_ligand_chemistry

        mol, traj, idx = self._ligand("c1ccccc1O")
        path = tmp_path / "given.sdf"
        Chem.SDWriter(str(path)).write(mol)

        got = resolve_ligand_chemistry(traj, "LIG", idx, supplied=path,
                                       allow_fetch=False)
        assert got.source == "supplied"
        assert not got.is_perceived

    def test_a_run_of_our_own_is_preferred_over_perceiving(self, tmp_path) -> None:
        from rdkit import Chem

        from fastmdxplora.analysis.ligand_chemistry import resolve_ligand_chemistry

        mol, traj, idx = self._ligand("c1ccccc1O")
        ligands = tmp_path / "setup" / "ligands"
        ligands.mkdir(parents=True)
        Chem.SDWriter(str(ligands / "LIG.sdf")).write(mol)

        got = resolve_ligand_chemistry(traj, "LIG", idx, run_dir=tmp_path,
                                       allow_fetch=False)
        assert got.source == "run"

    def test_a_definition_of_a_different_molecule_is_refused(
        self, tmp_path
    ) -> None:
        """Worse than none: it would put donors and rings on atoms that are
        not there."""
        from rdkit import Chem

        from fastmdxplora.analysis.ligand_chemistry import resolve_ligand_chemistry

        _mol, traj, idx = self._ligand("c1ccccc1O")
        other, _t, _i = self._ligand("CCO")
        path = tmp_path / "wrong.sdf"
        Chem.SDWriter(str(path)).write(other)

        got = resolve_ligand_chemistry(traj, "LIG", idx, supplied=path,
                                       allow_fetch=False)
        assert got.source == "perceived", "a mismatched definition was accepted"

    def test_perceiving_gets_ordinary_ligands_right(self) -> None:
        from fastmdxplora.analysis.ligand_chemistry import resolve_ligand_chemistry

        for smiles in ("c1ccccc1O", "NC(=N)c1ccccc1", "c1ccc2[nH]ccc2c1",
                       "Cn1cnc2c1c(=O)n(C)c(=O)n2C"):
            truth, traj, idx = self._ligand(smiles)
            got = resolve_ligand_chemistry(traj, "LIG", idx, allow_fetch=False)
            assert self._double_bonds(got.mol) == self._double_bonds(truth), smiles
            aromatic = sum(1 for a in got.mol.GetAtoms() if a.GetIsAromatic())
            expected = sum(1 for a in truth.GetAtoms() if a.GetIsAromatic())
            assert aromatic == expected, smiles

    def test_a_charge_that_only_one_value_explains_is_found(self) -> None:
        """Acetate is -1, and only -1 balances. Assuming neutral would give it
        no double bond at all -- and a carboxylate is what forms salt bridges,
        so this is where guessing wrong costs most."""
        from fastmdxplora.analysis.ligand_chemistry import resolve_ligand_chemistry

        truth, traj, idx = self._ligand("CC(=O)[O-]")
        got = resolve_ligand_chemistry(traj, "LIG", idx, allow_fetch=False)
        assert self._double_bonds(got.mol) == self._double_bonds(truth)
        assert not got.charge_was_ambiguous
        assert "-1" in got.detail

    def test_an_ambiguous_charge_is_reported_as_ambiguous(self) -> None:
        """Guanidinium is +1, and both -1 and +1 give a molecule that
        sanitises. Taking the first would make it an anion, the opposite of
        the charge that forms its salt bridges."""
        from fastmdxplora.analysis.ligand_chemistry import resolve_ligand_chemistry

        _truth, traj, idx = self._ligand("NC(=[NH2+])N")
        got = resolve_ligand_chemistry(traj, "LIG", idx, allow_fetch=False)
        assert got.charge_was_ambiguous
        assert "would also balance" in got.detail
        assert got.as_record()["charge_was_ambiguous"] is True

    def test_stating_the_charge_settles_it(self) -> None:
        from fastmdxplora.analysis.ligand_chemistry import resolve_ligand_chemistry

        truth, traj, idx = self._ligand("NC(=[NH2+])N")
        got = resolve_ligand_chemistry(traj, "LIG", idx, net_charge=1,
                                       allow_fetch=False)
        assert self._double_bonds(got.mol) == self._double_bonds(truth)
        assert not got.charge_was_ambiguous

    def test_the_results_carry_how_they_were_arrived_at(self) -> None:
        """So a reader knows whether an interaction rests on chemistry that
        was resolved or chemistry that was inferred."""
        from fastmdxplora.analysis.ligand_chemistry import resolve_ligand_chemistry

        _truth, traj, idx = self._ligand("c1ccccc1O")
        record = resolve_ligand_chemistry(
            traj, "LIG", idx, allow_fetch=False).as_record()
        assert record["source"] == "perceived"
        assert "guess" in record["confidence"]
        assert record["bond_orders_perceived"] is True
