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
        assert refusal_of(caught.value).details["names"] == ["AAA"]

    def test_a_name_is_compared_as_preparation_reads_it(self, tmp_path) -> None:
        # AAAX.sdf is built as AAA: its stem is cut to three letters.
        one = _write(tmp_path / "AAAX.sdf", SDF.format(title="AAAX"))
        many = _write(tmp_path / "set.sdf", SDF.format(title="AAA") + SDF.format(title="BBB"))
        with pytest.raises(LigandError) as caught:
            expand_ligand_files([str(one), str(many)], tmp_path / "ligands")
        assert refusal_of(caught.value).code == "setup.ligand.names_collide"

    def test_a_standard_residue_name_is_not_taken_as_a_title(self, tmp_path) -> None:
        # A ligand called ALA, HOH or NA would be selected with the protein,
        # the water or the ions.
        src = _write(tmp_path / "set.sdf", SDF.format(title="ALA") + SDF.format(title="hoh")
                     + SDF.format(title="NA") + SDF.format(title="BEN"))
        out = expand_ligand_files(str(src), tmp_path / "ligands")
        assert [Path(f).stem for f in out] == ["L01", "L02", "L03", "BEN"]

    def test_two_files_of_several_get_different_names(self, tmp_path) -> None:
        a = _write(tmp_path / "a.sdf", SDF.format(title="") * 2)
        b = _write(tmp_path / "b.sdf", SDF.format(title="") * 2)
        out = expand_ligand_files([str(a), str(b)], tmp_path / "ligands")
        assert [Path(f).stem for f in out] == ["L01", "L02", "L03", "L04"]
        assert _resolve_ligand_names(out, None) == ["L01", "L02", "L03", "L04"]

    def test_a_numbered_name_in_use_is_skipped(self, tmp_path) -> None:
        one = _write(tmp_path / "L01.sdf", SDF.format(title="X"))
        many = _write(tmp_path / "set.sdf", SDF.format(title="") + SDF.format(title="L03"))
        out = expand_ligand_files([str(one), str(many)], tmp_path / "ligands")
        assert [Path(f).stem for f in out] == ["L01", "L02", "L03"]

    def test_running_out_of_names_is_refused(self, tmp_path, monkeypatch) -> None:
        from fastmdxplora.setup import ligand as ligand_module

        monkeypatch.setattr(ligand_module, "_generated_names", lambda: iter(["L01"]))
        src = _write(tmp_path / "set.sdf", SDF.format(title="") * 2)
        with pytest.raises(LigandError) as caught:
            expand_ligand_files(str(src), tmp_path / "ligands")
        assert refusal_of(caught.value).code == "setup.ligand.too_many_molecules"

    def test_a_numbered_name_is_never_a_standard_residue(self, tmp_path, monkeypatch) -> None:
        from fastmdxplora.setup import ligand as ligand_module

        monkeypatch.setattr(ligand_module, "_generated_names",
                            lambda: iter(["LEU", "L01", "LYS", "L02"]))
        src = _write(tmp_path / "set.sdf", SDF.format(title="") * 2)
        out = expand_ligand_files(str(src), tmp_path / "ligands")
        assert [Path(f).stem for f in out] == ["L01", "L02"]

    def test_a_title_selected_as_protein_or_unselectable_is_not_taken(self, tmp_path) -> None:
        # HIP is a protonated histidine to MDTraj's `protein`; `resname 07`
        # is not a selection it can read.
        src = _write(tmp_path / "set.sdf", SDF.format(title="HIP") + SDF.format(title="07")
                     + SDF.format(title="B12"))
        out = expand_ligand_files(str(src), tmp_path / "ligands")
        assert [Path(f).stem for f in out] == ["L01", "L02", "B12"]

    def test_names_in_use_are_said_together(self, tmp_path) -> None:
        a = _write(tmp_path / "AAA.sdf", SDF.format(title="AAA"))
        b = _write(tmp_path / "BBB.sdf", SDF.format(title="BBB"))
        many = _write(tmp_path / "set.sdf", SDF.format(title="AAA") + SDF.format(title="BBB"))
        with pytest.raises(LigandError, match="are called AAA, BBB, names") as caught:
            expand_ligand_files([str(a), str(b), str(many)], tmp_path / "ligands")
        assert refusal_of(caught.value).details["names"] == ["AAA", "BBB"]

    def test_a_record_keeps_the_blank_line_that_ends_its_data(self, tmp_path) -> None:
        # The blank line before `$$$$` ends the last data item; without it a
        # reader takes `$$$$` into the item's value.
        first = SDF.format(title="AAA").replace("$$$$\n", "> <NOTE>\nbar\n\n$$$$  \n")
        src = _write(tmp_path / "set.sdf", first + SDF.format(title="BBB").rstrip("$\n"))
        out = expand_ligand_files(str(src), tmp_path / "ligands")
        assert Path(out[0]).read_text() == first
        assert Path(out[1]).read_text().endswith("M  END\n$$$$\n")

    def test_names_given_by_the_caller_make_no_clash(self, tmp_path) -> None:
        one = _write(tmp_path / "AAA.sdf", SDF.format(title="AAA"))
        many = _write(tmp_path / "set.sdf", SDF.format(title="AAA") + SDF.format(title="BBB"))
        out = expand_ligand_files([str(one), str(many)], tmp_path / "ligands", named=True)
        assert [Path(f).stem for f in out] == ["AAA", "L01", "BBB"]

    def test_a_shared_title_never_used_makes_no_clash(self, tmp_path) -> None:
        one = _write(tmp_path / "AAA.sdf", SDF.format(title="AAA"))
        many = _write(tmp_path / "set.sdf", SDF.format(title="AAA") * 2)
        out = expand_ligand_files([str(one), str(many)], tmp_path / "ligands")
        assert [Path(f).stem for f in out] == ["AAA", "L01", "L02"]

    def test_a_title_in_another_encoding_is_told_readably(self, tmp_path) -> None:
        src = tmp_path / "set.sdf"
        src.write_bytes((SDF.format(title="Caf\xe9") + SDF.format(title="BBB")).encode("latin-1"))
        told: list[str] = []
        expand_ligand_files(str(src), tmp_path / "ligands", told=told)
        assert "'Caf\ufffd' as L01" in told[0]

    def test_the_numbered_names_go_past_ninety_nine(self) -> None:
        from fastmdxplora.setup.ligand import _NAME, _generated_names

        names = list(_generated_names())
        assert names[:2] == ["L01", "L02"] and names[98:100] == ["L99", "LA0"]
        assert len(set(names)) == len(names) and all(_NAME.fullmatch(n) for n in names)

    def test_a_record_is_written_byte_for_byte(self, tmp_path) -> None:
        # A Latin-1 name in a comment line and Windows line endings both
        # reach the toolkit as they were written.
        body = SDF.replace("\n", "\r\n")
        raw = (body.format(title="AAA").replace("\r\n\r\n", "\r\nCaf\xe9\r\n", 1)
               + body.format(title="BBB")).encode("latin-1")
        src = tmp_path / "set.sdf"
        src.write_bytes(raw)
        out = expand_ligand_files(str(src), tmp_path / "ligands")
        first = Path(out[0]).read_bytes()
        assert b"Caf\xe9" in first
        assert b"\r\n" in first and b"\n" not in first.replace(b"\r\n", b"")
        assert first.rstrip(b"\r\n").endswith(b"$$$$")
        assert first.count(b"$$$$") == 1

    def test_the_division_is_told(self, tmp_path) -> None:
        src = _write(tmp_path / "set.sdf", SDF.format(title="benzamidine") + SDF.format(title="ATP"))
        told: list[str] = []
        expand_ligand_files(str(src), tmp_path / "ligands", told=told)
        assert told == ["The ligand file set.sdf holds 2 molecules, each prepared as "
                        "its own ligand: 'benzamidine' as L01, ATP."]


class TestSetupHandsPreparationOneMoleculePerFile:
    """The division happens in setup, before preparation sees the ligands."""

    def _run(self, tmp_path, monkeypatch, structure_extra: str = "", **options):
        from unittest.mock import MagicMock, patch

        from fastmdxplora.setup import ligand as ligand_module
        from fastmdxplora.setup import pdbfix as pdbfix_module
        from fastmdxplora.setup import pipeline
        from fastmdxplora.setup import prepare as prepare_module
        from tests.test_a_ligand_needs_a_charge_provider import TRIPEPTIDE

        structure = tmp_path / "complex.pdb"
        structure.write_text(TRIPEPTIDE + structure_extra, encoding="utf-8")
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
        # The default name is not a name given to both of them; each is
        # named by its file.
        assert kwargs["ligand_name"] == ["AAA", "BBB"]

    def test_a_name_given_by_the_caller_is_kept(self, tmp_path, monkeypatch) -> None:
        src = _write(tmp_path / "set.sdf", SDF.format(title="AAA") + SDF.format(title="BBB"))
        solvate = self._run(tmp_path, monkeypatch, ligand=str(src),
                            ligand_name=["XX1", "XX2"])
        assert solvate.call_args.kwargs["ligand_name"] == ["XX1", "XX2"]

    def test_one_charge_for_several_molecules_is_refused(self, tmp_path, monkeypatch) -> None:
        from fastmdxplora.refusals import StudyError

        src = _write(tmp_path / "set.sdf", SDF.format(title="AAA") + SDF.format(title="BBB"))
        with pytest.raises(StudyError, match="each molecule's own formal charges") as caught:
            self._run(tmp_path, monkeypatch, ligand=str(src), ligand_net_charge=1)
        refusal = refusal_of(caught.value)
        assert refusal.code == "config.option.conflicting"
        assert refusal.details["options"] == ["ligand_net_charge", "ligand"]

    def test_a_charge_per_molecule_is_accepted(self, tmp_path, monkeypatch) -> None:
        src = _write(tmp_path / "set.sdf", SDF.format(title="AAA") + SDF.format(title="BBB"))
        solvate = self._run(tmp_path, monkeypatch, ligand=str(src),
                            ligand_net_charge=[1, 0])
        assert solvate.call_args.kwargs["ligand_net_charge"] == [1, 0]

    def test_two_files_need_no_name_either(self, tmp_path, monkeypatch) -> None:
        a = _write(tmp_path / "AAA.sdf", SDF.format(title="AAA"))
        b = _write(tmp_path / "BBB.sdf", SDF.format(title="BBB"))
        solvate = self._run(tmp_path, monkeypatch, ligand=[str(a), str(b)])
        assert solvate.call_args.kwargs["ligand_name"] == ["AAA", "BBB"]

    def test_the_manifest_records_the_names_used(self, tmp_path, monkeypatch) -> None:
        import json

        src = _write(tmp_path / "set.sdf", SDF.format(title="AAA") + SDF.format(title=""))
        self._run(tmp_path, monkeypatch, ligand=str(src))
        manifests = list((tmp_path / "setup").rglob("*.json"))
        recorded = [json.loads(m.read_text()) for m in manifests]
        names = [r["resolved_forcefield"]["ligand"]["name"] for r in recorded
                 if "ligand" in (r.get("resolved_forcefield") or {})]
        assert names == [["AAA", "L01"]]

    def test_the_default_name_read_back_is_not_a_name_given(self, tmp_path, monkeypatch) -> None:
        # A study's resolved config records the default; read back, it comes
        # under both spellings.
        src = _write(tmp_path / "set.sdf", SDF.format(title="AAA") + SDF.format(title="BBB"))
        solvate = self._run(tmp_path, monkeypatch, ligand=str(src),
                            ligand_name="LIG", ligand_resname="LIG")
        assert solvate.call_args.kwargs["ligand_name"] == ["AAA", "BBB"]

    @staticmethod
    def _parameters(tmp_path) -> dict:
        import json

        return json.loads((tmp_path / "setup" / "setup_parameters.json").read_text())

    def test_the_setup_record_names_the_ligands(self, tmp_path, monkeypatch) -> None:
        # What the methods and the repair read: the names, not the files.
        src = _write(tmp_path / "set.sdf", SDF.format(title="AAA") + SDF.format(title=""))
        self._run(tmp_path, monkeypatch, ligand=str(src))
        parameters = self._parameters(tmp_path)["parameters"]
        assert parameters["ligand_name"] == ["AAA", "L01"]
        assert parameters["_reinstated_heterogens"] == ["AAA", "L01"]

    def test_a_molecule_not_named_for_the_structure_s_ligand_is_said(
            self, tmp_path, monkeypatch) -> None:
        ben = ("HETATM  100  C1  BEN B 401       1.000   1.000   1.000"
               "  1.00  0.00           C\n")
        src = _write(tmp_path / "set.sdf",
                     SDF.format(title="benzamidine") + SDF.format(title="ATP"))
        self._run(tmp_path, monkeypatch, structure_extra=ben, ligand=str(src))
        notes = self._parameters(tmp_path)["notes"]
        said = [n for n in notes if "take their poses from the file" in n]
        assert said and "holds BEN" in said[0] and "L01, ATP take" in said[0]

    @staticmethod
    def _hetatm(name: str) -> str:
        return (f"HETATM  100  C1  {name:>3} B 401       1.000   1.000   1.000"
                "  1.00  0.00           C\n")

    def test_a_cofactor_a_title_cannot_name_is_said_with_its_remedy(
            self, tmp_path, monkeypatch) -> None:
        # SAH is a cofactor MDTraj selects as protein, and 7VF starts with a
        # digit: neither is taken from a title, so the molecule is L01 and
        # the way to place it is a file of its own.
        for code in ("SAH", "7VF"):
            (tmp_path / code).mkdir()
            src = _write(tmp_path / "set.sdf", SDF.format(title=code) + SDF.format(title="XYZ"))
            self._run(tmp_path / code, monkeypatch, structure_extra=self._hetatm(code),
                      ligand=str(src))
            notes = self._parameters(tmp_path / code)["notes"]
            said = [n for n in notes if "from the file" in n]
            assert said and f"holds {code}" in said[0]
            assert f"A title cannot name {code}" in said[0]
            assert f"({code}.sdf)" in said[0] and "as its title" not in said[0]

    def test_one_molecule_takes_its_pose(self, tmp_path, monkeypatch) -> None:
        src = _write(tmp_path / "set.sdf", SDF.format(title="benzamidine") + SDF.format(title="BEN"))
        self._run(tmp_path, monkeypatch,
                  structure_extra=self._hetatm("BEN") + self._hetatm("ATP"), ligand=str(src))
        said = [n for n in self._parameters(tmp_path)["notes"] if "from the file" in n]
        assert said and "holds ATP" in said[0] and "L01 takes its pose" in said[0]
        assert "the structure's ATP, give it" in said[0]

    def test_a_ligand_supplied_as_its_own_file_is_not_warned_for(
            self, tmp_path, monkeypatch) -> None:
        ben = _write(tmp_path / "BEN.sdf", SDF.format(title="BEN"))
        src = _write(tmp_path / "set.sdf", SDF.format(title="") * 2)
        self._run(tmp_path, monkeypatch, structure_extra=self._hetatm("BEN"),
                  ligand=[str(ben), str(src)])
        assert not [n for n in self._parameters(tmp_path)["notes"] if "from the file" in n]

    def test_names_given_by_the_caller_are_the_names_checked(
            self, tmp_path, monkeypatch) -> None:
        src = _write(tmp_path / "set.sdf", SDF.format(title="") * 2)
        self._run(tmp_path, monkeypatch, structure_extra=self._hetatm("BEN"),
                  ligand=str(src), ligand_name=["BEN", "XYZ"])
        assert not [n for n in self._parameters(tmp_path)["notes"] if "from the file" in n]

    def test_poses_asked_for_from_the_file_are_not_warned_for(
            self, tmp_path, monkeypatch) -> None:
        src = _write(tmp_path / "set.sdf", SDF.format(title="") * 2)
        self._run(tmp_path, monkeypatch, structure_extra=self._hetatm("BEN"),
                  ligand=str(src), ligand_pose="FILE")
        assert not [n for n in self._parameters(tmp_path)["notes"] if "from the file" in n]

    def test_a_molecule_named_for_the_structure_s_ligand_is_not_warned(
            self, tmp_path, monkeypatch) -> None:
        ben = ("HETATM  100  C1  BEN B 401       1.000   1.000   1.000"
               "  1.00  0.00           C\n")
        src = _write(tmp_path / "set.sdf", SDF.format(title="BEN") + SDF.format(title="ATP"))
        self._run(tmp_path, monkeypatch, structure_extra=ben, ligand=str(src))
        notes = self._parameters(tmp_path)["notes"]
        assert not [n for n in notes if "take their poses from the file" in n]


class TestSeveralLigandsAreNamedDownstream:
    @staticmethod
    def _manifest(root: Path, name) -> None:
        import json

        (root / "setup").mkdir(parents=True)
        (root / "setup" / "setup_parameters.json").write_text(json.dumps(
            {"resolved_forcefield": {"ligand": {"name": name}}}), encoding="utf-8")

    def test_analysis_selects_every_ligand(self, tmp_path) -> None:
        from fastmdxplora.analysis.analyze import _detect_ligand_resname

        self._manifest(tmp_path, ["AAA", "L01"])
        assert _detect_ligand_resname(tmp_path) == "AAA L01"

    def test_the_cross_tool_check_selects_every_ligand(self, tmp_path) -> None:
        from fastmdxplora.validation.cross_tool import (
            _configured_ligand_resname, topology_has_resname)

        self._manifest(tmp_path, ["AAA", "L01"])
        assert _configured_ligand_resname(tmp_path) == "AAA L01"
        top = tmp_path / "top.pdb"
        top.write_text("HETATM    1  C1  L01 A   1       0.000   0.000   0.000"
                       "  1.00  0.00           C\n", encoding="utf-8")
        assert topology_has_resname(top, "AAA L01")
        assert not topology_has_resname(top, "BBB")
