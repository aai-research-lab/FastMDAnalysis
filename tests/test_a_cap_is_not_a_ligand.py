"""Capping groups belong to the polymer, so nothing else claims them.

Alanine dipeptide reported one protein residue, ten protein atoms and two
ligands. The molecule is three residues of twenty-two atoms and has no ligand
at all: ACE and NME are the caps that give its single alanine proper backbone
neighbours.

Worse than a miscount. A cap surfaced as a ligand invites the protein-ligand
analyses -- contacts, hydrogen bonds, interaction fingerprints -- between a
residue and its own backbone terminus, and a reader of that report would be
looking at a study of nothing.

The setup phase already knew: heterogen removal keeps these because they are
the molecule. What did not know was everything downstream.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from fastmdxplora.gui.ligand_detection import (
    AMINO_ACID_RESNAMES, CAPPING_RESNAMES, POLYMER_RESNAMES, detect_ligands)

_CAPPED_ALANINE = [
    ("HH31", "ACE", 1), ("CH3", "ACE", 1), ("HH32", "ACE", 1),
    ("HH33", "ACE", 1), ("C", "ACE", 1), ("O", "ACE", 1),
    ("N", "ALA", 2), ("H", "ALA", 2), ("CA", "ALA", 2), ("HA", "ALA", 2),
    ("CB", "ALA", 2), ("HB1", "ALA", 2), ("HB2", "ALA", 2),
    ("HB3", "ALA", 2), ("C", "ALA", 2), ("O", "ALA", 2),
    ("N", "NME", 3), ("H", "NME", 3), ("CH3", "NME", 3),
    ("HH31", "NME", 3), ("HH32", "NME", 3), ("HH33", "NME", 3),
]


def _write_capped_alanine(directory):
    """Three residues, twenty-two atoms, no ligand.

    Written here rather than read from a file beside the repository: the
    first version skipped wherever that file was absent, which was
    everywhere except the machine it was written on.
    """
    lines = []
    for index, (name, residue, sequence) in enumerate(_CAPPED_ALANINE, 1):
        atom = name.ljust(4) if len(name) == 4 else (" " + name).ljust(4)
        lines.append(
            f"ATOM  {index:5d} {atom} {residue:>3s} A{sequence:4d}    "
            f"{index:8.3f}{index:8.3f}{index:8.3f}  1.00  0.00"
            f"          {name[0]:>2s}")
    lines.append("END")
    path = directory / "capped-alanine.pdb"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


class TestTheCapsAreCountedAsPolymer:
    def test_the_usual_ones_are_listed(self) -> None:
        assert {"ACE", "NME", "NHE"} <= CAPPING_RESNAMES

    def test_the_polymer_is_the_two_together(self) -> None:
        assert POLYMER_RESNAMES == AMINO_ACID_RESNAMES | CAPPING_RESNAMES

    def test_a_cap_is_not_offered_as_a_ligand(self) -> None:
        keys = [("A", "ACE", "1"), ("A", "ALA", "2"), ("A", "NME", "3")]
        assert detect_ligands(keys)["resnames"] == []

    def test_the_counts_describe_the_molecule(self, tmp_path) -> None:
        from fastmdxplora.gui.structure_info import count_structure

        got = count_structure(_write_capped_alanine(tmp_path))
        assert got["protein_residues"] == 3
        assert got["protein_atoms"] == 22
        assert got["ligand_resnames"] == []


class TestARealLigandStillSurfaces:
    """The point is a short known list, not a licence to swallow everything
    that is not an amino acid."""

    def test_a_bound_small_molecule_is_a_ligand(self) -> None:
        keys = [("A", "ALA", "1"), ("A", "BNZ", "2")]
        assert detect_ligands(keys)["resnames"] == ["BNZ"]

    def test_a_cap_beside_a_ligand_does_not_hide_it(self) -> None:
        keys = [("A", "ACE", "1"), ("A", "ALA", "2"), ("A", "NME", "3"),
                ("A", "BNZ", "4")]
        assert detect_ligands(keys)["resnames"] == ["BNZ"]

    def test_an_explicit_name_still_wins(self) -> None:
        """Somebody who says a cap is their ligand is not overruled."""
        keys = [("A", "ALA", "1"), ("A", "ACE", "2")]
        assert detect_ligands(keys, explicit={"ACE"})["resnames"] == ["ACE"]


class TestTheTwoConsumersAgree:
    def test_counting_and_detecting_read_the_same_list(self) -> None:
        """Two ideas of what the polymer is would drift, and the drift would
        show up as a residue counted out of the protein and in as a ligand --
        which is the bug this file is about."""
        import fastmdxplora.gui.structure_info as info

        source = Path(info.__file__).read_text(encoding="utf-8")
        assert "POLYMER_RESNAMES" in source
        assert "in AMINO_ACID_RESNAMES" not in source
