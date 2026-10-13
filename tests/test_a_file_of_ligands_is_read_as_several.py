"""A file holding several molecules is several ligands.

A list of N ligand files has always been N ligands. One file holding N
molecules is the same study written differently, and used to be refused.
"""

from pathlib import Path

import pytest

from fastmdxplora.refusals import refusal_of
from fastmdxplora.setup.ligand import LigandError, expand_ligand_files
from fastmdxplora.setup.prepare import _resolve_ligand_names

SDF = "{title}\n\n\n  0  0  0  0  0  0  0  0  0  0999 V2000\nM  END\n$$$$\n"
MOL2 = "@<TRIPOS>MOLECULE\n{title}\n 0 0 0 0 0\nSMALL\nNO_CHARGES\n\n"


def _write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


class TestAFileOfSeveralMolecules:
    def test_an_sdf_is_divided_one_file_per_molecule(self, tmp_path) -> None:
        src = _write(tmp_path / "set.sdf", SDF.format(title="BEN") + SDF.format(title="ATP"))
        out = expand_ligand_files(str(src), tmp_path / "ligands")
        assert [Path(f).name for f in out] == ["BEN.sdf", "ATP.sdf"]
        assert all(Path(f).read_text().count("$$$$") == 1 for f in out)

    def test_a_mol2_is_divided_the_same_way(self, tmp_path) -> None:
        src = _write(tmp_path / "set.mol2", MOL2.format(title="AAA") + MOL2.format(title="BBB"))
        out = expand_ligand_files([str(src)], tmp_path / "ligands")
        assert [Path(f).name for f in out] == ["AAA.mol2", "BBB.mol2"]

    def test_unusable_titles_become_numbered_names(self, tmp_path) -> None:
        src = _write(tmp_path / "set.sdf",
                     SDF.format(title="benzamidine hydrochloride") + SDF.format(title=""))
        out = expand_ligand_files(str(src), tmp_path / "ligands")
        assert [Path(f).stem for f in out] == ["L01", "L02"]

    def test_a_title_shared_by_two_molecules_is_not_trusted(self, tmp_path) -> None:
        src = _write(tmp_path / "set.sdf", SDF.format(title="LIG") + SDF.format(title="LIG"))
        out = expand_ligand_files(str(src), tmp_path / "ligands")
        assert [Path(f).stem for f in out] == ["L01", "L02"]

    def test_the_new_names_are_distinct_residue_names(self, tmp_path) -> None:
        src = _write(tmp_path / "set.sdf", SDF.format(title="") * 3)
        out = expand_ligand_files(str(src), tmp_path / "ligands")
        assert _resolve_ligand_names(out, None) == ["L01", "L02", "L03"]

    def test_a_single_molecule_file_is_returned_untouched(self, tmp_path) -> None:
        src = _write(tmp_path / "one.sdf", SDF.format(title="BEN"))
        given = str(src)
        assert expand_ligand_files(given, tmp_path / "ligands") is given
        assert not (tmp_path / "ligands").exists()

    def test_a_missing_file_is_left_for_the_loader_to_describe(self, tmp_path) -> None:
        given = [str(tmp_path / "absent.sdf")]
        assert expand_ligand_files(given, tmp_path / "ligands") is given

    def test_a_file_beside_another_keeps_its_place(self, tmp_path) -> None:
        one = _write(tmp_path / "ATP.sdf", SDF.format(title="ATP"))
        many = _write(tmp_path / "set.sdf", SDF.format(title="AAA") + SDF.format(title="BBB"))
        out = expand_ligand_files([str(one), str(many)], tmp_path / "ligands")
        assert [Path(f).stem for f in out] == ["ATP", "AAA", "BBB"]

    def test_a_name_already_in_use_is_refused(self, tmp_path) -> None:
        one = _write(tmp_path / "AAA.sdf", SDF.format(title="AAA"))
        many = _write(tmp_path / "set.sdf", SDF.format(title="AAA") + SDF.format(title="BBB"))
        with pytest.raises(LigandError) as caught:
            expand_ligand_files([str(one), str(many)], tmp_path / "ligands")
        assert refusal_of(caught.value).code == "setup.ligand.names_collide"


class TestSetupHandsPreparationOneMoleculePerFile:
    """The division happens in setup, before preparation sees the ligands."""

    def _run(self, tmp_path, monkeypatch, **options):
        from unittest.mock import MagicMock, patch

        from fastmdxplora.setup import ligand as ligand_module
        from fastmdxplora.setup import pdbfix as pdbfix_module
        from fastmdxplora.setup import pipeline
        from fastmdxplora.setup import prepare as prepare_module
        from tests.test_a_ligand_needs_a_charge_provider import TRIPEPTIDE

        structure = tmp_path / "complex.pdb"
        structure.write_text(TRIPEPTIDE, encoding="utf-8")
        monkeypatch.setattr(ligand_module, "am1bcc_provider", lambda: "AmberTools")
        orchestrator = MagicMock()
        orchestrator.system = str(structure)
        orchestrator._presenter = None
        setup_dir = tmp_path / "setup"
        setup_dir.mkdir()
        with patch.object(pdbfix_module, "fix_pdb_with_pdbfixer"), \
                patch.object(prepare_module, "prepare_system",
                             return_value={}) as solvate:
            pipeline.run(orchestrator=orchestrator, output_dir=setup_dir, **options)
        return solvate

    def test_preparation_is_given_the_molecules_separately(self, tmp_path, monkeypatch) -> None:
        src = _write(tmp_path / "set.sdf", SDF.format(title="AAA") + SDF.format(title="BBB"))
        solvate = self._run(tmp_path, monkeypatch, ligand=str(src))
        kwargs = solvate.call_args.kwargs
        assert [Path(f).name for f in kwargs["ligand"]] == ["AAA.sdf", "BBB.sdf"]
        # The default name is not a name given to both of them.
        assert kwargs["ligand_name"] is None

    def test_a_name_given_by_the_caller_is_kept(self, tmp_path, monkeypatch) -> None:
        src = _write(tmp_path / "set.sdf", SDF.format(title="AAA") + SDF.format(title="BBB"))
        solvate = self._run(tmp_path, monkeypatch, ligand=str(src),
                            ligand_name=["XX1", "XX2"])
        assert solvate.call_args.kwargs["ligand_name"] == ["XX1", "XX2"]

    def test_one_charge_reaches_every_molecule_as_it_does_for_separate_files(
            self, tmp_path, monkeypatch) -> None:
        """A stated charge is checked against each molecule's own formal
        charges and never substituted for them, so one number is as safe for
        an expanded file as it is for a list of files."""
        src = _write(tmp_path / "set.sdf", SDF.format(title="AAA") + SDF.format(title="BBB"))
        solvate = self._run(tmp_path, monkeypatch, ligand=str(src), ligand_net_charge=0)
        assert solvate.call_args.kwargs["ligand_net_charge"] == 0

    def test_a_charge_per_molecule_is_accepted(self, tmp_path, monkeypatch) -> None:
        src = _write(tmp_path / "set.sdf", SDF.format(title="AAA") + SDF.format(title="BBB"))
        solvate = self._run(tmp_path, monkeypatch, ligand=str(src),
                            ligand_net_charge=[1, 0])
        assert solvate.call_args.kwargs["ligand_net_charge"] == [1, 0]

    def test_two_files_need_no_name_either(self, tmp_path, monkeypatch) -> None:
        a = _write(tmp_path / "AAA.sdf", SDF.format(title="AAA"))
        b = _write(tmp_path / "BBB.sdf", SDF.format(title="BBB"))
        solvate = self._run(tmp_path, monkeypatch, ligand=[str(a), str(b)])
        assert solvate.call_args.kwargs["ligand_name"] is None
