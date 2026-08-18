"""The weakest chemistry should not be the only kind with no record.

Chemistry read from a file is already on disk. Chemistry inferred from
coordinates is not, and it is the route the software itself labels "bond
orders are a guess". A benchmark run that perceived its ligand left no
trace of what it had perceived, so every later analysis inferred it again
from scratch -- and if the inference ever differed, nothing would show it.
"""

from __future__ import annotations

import pytest


def _perceived(resname="BEN"):
    from rdkit import Chem

    from fastmdxplora.analysis.ligand_chemistry import ResolvedChemistry

    mol = Chem.AddHs(Chem.MolFromSmiles("NC(=[NH2+])c1ccccc1"))
    return ResolvedChemistry(
        mol=mol, source="perceived",
        detail="bond orders inferred from the coordinates",
        resname=resname, n_atoms=mol.GetNumAtoms())


def _from_file(resname="BEN"):
    from rdkit import Chem

    from fastmdxplora.analysis.ligand_chemistry import ResolvedChemistry

    mol = Chem.AddHs(Chem.MolFromSmiles("NC(=[NH2+])c1ccccc1"))
    return ResolvedChemistry(
        mol=mol, source="supplied", detail="from an SDF",
        resname=resname, n_atoms=mol.GetNumAtoms())


class TestWhatGetsDeposited:
    def test_a_perceived_molecule_is_written(self, tmp_path):
        from fastmdxplora.analysis.ligand_chemistry import (
            deposit_perceived_chemistry,
        )

        pytest.importorskip("rdkit", reason="requires the [ligand] extra")
        written = deposit_perceived_chemistry(_perceived(), tmp_path)

        assert written is not None and written.is_file()
        assert written.name == "BEN.sdf"
        assert written.parent.name == "ligands"

    def test_it_lands_where_resolution_looks(self, tmp_path):
        """The search is `**/ligands/<resname>*.sdf`, so the next analysis
        finds it as a resolved file rather than inferring again."""
        pytest.importorskip("rdkit", reason="requires the [ligand] extra")
        from fastmdxplora.analysis.ligand_chemistry import (
            deposit_perceived_chemistry,
        )

        deposit_perceived_chemistry(_perceived(), tmp_path)
        assert list(tmp_path.glob("**/ligands/BEN*.sdf"))

    def test_chemistry_from_a_file_is_not_rewritten(self, tmp_path):
        """It is already on disk; a second copy is a second thing to keep
        in step."""
        pytest.importorskip("rdkit", reason="requires the [ligand] extra")
        from fastmdxplora.analysis.ligand_chemistry import (
            deposit_perceived_chemistry,
        )

        assert deposit_perceived_chemistry(_from_file(), tmp_path) is None

    def test_an_existing_deposit_is_left_alone(self, tmp_path):
        pytest.importorskip("rdkit", reason="requires the [ligand] extra")
        from fastmdxplora.analysis.ligand_chemistry import (
            deposit_perceived_chemistry,
        )

        first = deposit_perceived_chemistry(_perceived(), tmp_path)
        first.write_text("MARKER", encoding="utf-8")
        again = deposit_perceived_chemistry(_perceived(), tmp_path)

        assert again == first
        assert first.read_text(encoding="utf-8") == "MARKER", (
            "a deposit must not overwrite what a run already recorded")

    def test_no_run_directory_is_not_an_error(self):
        """An analysis run outside a project still has its chemistry; only
        the record is lost, and losing a run to save a record is the wrong
        trade."""
        pytest.importorskip("rdkit", reason="requires the [ligand] extra")
        from fastmdxplora.analysis.ligand_chemistry import (
            deposit_perceived_chemistry,
        )

        assert deposit_perceived_chemistry(_perceived(), None) is None

    def test_the_deposit_does_not_claim_certainty(self, tmp_path):
        """The record beside it still says perceived.

        Depositing what was inferred is a statement of what was used, not
        a promise that it was right.
        """
        pytest.importorskip("rdkit", reason="requires the [ligand] extra")
        from fastmdxplora.analysis.ligand_chemistry import (
            deposit_perceived_chemistry,
        )

        chemistry = _perceived()
        deposit_perceived_chemistry(chemistry, tmp_path)
        assert chemistry.as_record()["source"] == "perceived"
