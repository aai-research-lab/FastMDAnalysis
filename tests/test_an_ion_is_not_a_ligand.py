"""A lone ion has no chemistry to fetch.

An SDF carries bond orders, formal charges and aromaticity -- what a PDB
cannot express and a force field needs for a small molecule. A monatomic ion
has none of them, and the protein force field carries its parameters
directly.

Asked for anyway, the refusal told somebody to fetch `ZN_ideal.sdf`: a file
describing bonds that do not exist. 4INS is insulin with two structural zincs
and could not be prepared at all -- not because anything was wrong with it,
but because the pipeline asked a question with no answer.

Decided from the structure, not from a list of ion names: a list would need
every ion anybody simulates and would be wrong for whichever was missing.
"""

from __future__ import annotations

import logging
import tempfile
from pathlib import Path

import pytest

from fastmdxplora.setup.heterogens import Action, resolve
from fastmdxplora.setup.pipeline import _auto_ligands


def _structure(tmp_path: Path, *, ion: bool = True, ligand: bool = False,
               ion_name: str = "ZN") -> Path:
    """Histidines, optionally an ion on one, optionally an organic ligand."""
    lines, n = [], 0
    for chain, x in (("A", 10.0), ("B", 40.0)):
        n += 1
        lines.append(f"ATOM  {n:5d}  NE2 HIS {chain}  10    "
                     f"{x:8.3f}{10.0:8.3f}{10.0:8.3f}  1.00  0.00           N")
    if ion:
        n += 1
        lines.append(f"HETATM{n:5d} {ion_name:<4s} {ion_name:>3s} A 401    "
                     f"{10.2:8.3f}{10.0:8.3f}{10.0:8.3f}  1.00  0.00          "
                     f"{ion_name:>2s}")
    if ligand:
        for i, dx in enumerate((0.0, 1.4), 1):
            n += 1
            lines.append(f"HETATM{n:5d}  C{i}  BNZ A 500    "
                         f"{40.2 + dx:8.3f}{10.0:8.3f}{10.0:8.3f}"
                         f"  1.00  0.00           C")
    path = tmp_path / "system.pdb"
    path.write_text("\n".join(lines) + "\nEND\n", encoding="utf-8")
    return path


class TestTheStructureSaysWhichIsAnIon:
    def test_one_atom_is_monatomic(self, tmp_path: Path) -> None:
        decisions = resolve(_structure(tmp_path), keep_water=False)
        zinc = [d for d in decisions if d.resname == "ZN"]
        assert zinc and zinc[0].is_monatomic

    def test_several_atoms_are_not(self, tmp_path: Path) -> None:
        decisions = resolve(_structure(tmp_path, ion=False, ligand=True),
                            keep_water=False)
        benzene = [d for d in decisions if d.resname == "BNZ"]
        assert benzene and not benzene[0].is_monatomic

    @pytest.mark.parametrize("name", ["ZN", "MG", "CA", "FE", "NI"])
    def test_it_is_not_a_list_of_names(self, tmp_path: Path, name: str) -> None:
        """Any lone atom, whatever it is called. A list would be wrong for
        whichever ion somebody simulates that is not on it."""
        decisions = resolve(_structure(tmp_path, ion_name=name),
                            keep_water=False)
        found = [d for d in decisions if d.resname == name]
        assert found and found[0].is_monatomic


class TestAnIonIsKeptWithoutBeingFetched:
    def test_a_local_structure_with_only_ions_is_not_refused(
            self, tmp_path: Path) -> None:
        """The whole failure: 4INS could not be prepared at all."""
        assert _auto_ligands({"heterogens": "auto"},
                             _structure(tmp_path), tmp_path, None) == []

    def test_it_says_what_it_kept_and_why(self, tmp_path: Path, caplog) -> None:
        with caplog.at_level(logging.INFO):
            _auto_ligands({"heterogens": "auto"}, _structure(tmp_path),
                          tmp_path, None)
        assert "no chemistry to retrieve" in caplog.text
        assert "ZN" in caplog.text

    def test_the_plural_reads(self, tmp_path: Path, caplog) -> None:
        with caplog.at_level(logging.INFO):
            _auto_ligands({"heterogens": "auto"}, _structure(tmp_path),
                          tmp_path, None)
        assert "as an ion" in caplog.text


class TestAnOrganicLigandIsUnchanged:
    def test_it_still_refuses_without_an_entry_to_look_up(
            self, tmp_path: Path) -> None:
        with pytest.raises(ValueError, match="BNZ"):
            _auto_ligands({"heterogens": "auto"},
                          _structure(tmp_path, ligand=True), tmp_path, None)

    def test_the_ion_is_not_named_in_that_refusal(self, tmp_path: Path) -> None:
        """Naming it would send somebody after `ZN_ideal.sdf`, which is the
        fault this file is about."""
        with pytest.raises(ValueError) as caught:
            _auto_ligands({"heterogens": "auto"},
                          _structure(tmp_path, ligand=True), tmp_path, None)
        first = str(caught.value).split("\n")[0]
        assert "BNZ" in first and "ZN" not in first
