"""Point mutations, and the guard that makes them safe to write.

A mutation is stated against a numbering, and numbering travels badly:
the same protein appears in the literature with the construct's numbering,
the deposition's, and the mature sequence's. PDBFixer applies what it is
told, so a study asking for L99A against the wrong numbering would replace
whatever sits at 99 and every result afterwards would describe a protein
nobody meant to simulate.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from fastmdxplora.setup.pdbfix import (
    ONE_TO_THREE,
    fix_pdb_with_pdbfixer,
    parse_mutation,
)

TRIPEPTIDE = """\
ATOM      1  N   MET A   1       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  MET A   1       1.458   0.000   0.000  1.00  0.00           C
ATOM      3  C   MET A   1       2.009   1.420   0.000  1.00  0.00           C
ATOM      4  O   MET A   1       1.251   2.390   0.000  1.00  0.00           O
ATOM      5  CB  MET A   1       1.988  -0.773  -1.199  1.00  0.00           C
ATOM      6  N   LEU A   2       3.332   1.552   0.000  1.00  0.00           N
ATOM      7  CA  LEU A   2       3.999   2.845   0.000  1.00  0.00           C
ATOM      8  C   LEU A   2       5.510   2.664   0.000  1.00  0.00           C
ATOM      9  O   LEU A   2       6.020   1.545   0.000  1.00  0.00           O
ATOM     10  CB  LEU A   2       3.576   3.675  -1.213  1.00  0.00           C
ATOM     11  N   GLY A   3       6.235   3.777   0.000  1.00  0.00           N
ATOM     12  CA  GLY A   3       7.689   3.762   0.000  1.00  0.00           C
ATOM     13  C   GLY A   3       8.230   5.183   0.000  1.00  0.00           C
ATOM     14  O   GLY A   3       7.472   6.153   0.000  1.00  0.00           O
TER
END
"""


@pytest.fixture()
def structure(tmp_path) -> Path:
    path = tmp_path / "in.pdb"
    path.write_text(TRIPEPTIDE, encoding="utf-8")
    return path


class TestReadingAMutation:
    @pytest.mark.parametrize("text", ["L99A", "l99a", "LEU-99-ALA"])
    def test_both_conventions_read_the_same(self, text):
        assert parse_mutation(text) == ("LEU", 99, "ALA")

    @pytest.mark.parametrize(
        "text", ["99A", "LZ9A", "L99", "LEU-x-ALA", "LEU-99"])
    def test_what_is_not_a_mutation_is_refused(self, text):
        with pytest.raises(ValueError):
            parse_mutation(text)

    def test_a_letter_outside_the_twenty_is_refused(self):
        """Rather than guessed at.

        PDBFixer asks for three letters to avoid ambiguity, and this
        refuses for the same reason: a one-letter code it cannot place is
        not a residue it should invent a name for.
        """
        assert "X" not in ONE_TO_THREE
        with pytest.raises(ValueError, match="twenty amino acids"):
            parse_mutation("X99A")


class TestApplyingOne:
    """These need the engine; the parsing tests above deliberately do not.

    Reading a mutation and checking it against a structure is arithmetic
    and string handling, so it runs on every leg of the matrix. Applying
    one needs PDBFixer, which is in the `[md]` extra and is not installed
    on the floor leg. Skipped by name rather than by the import failing at
    call time, so the skip appears in the report: an analysis that reads
    as tested while its tests quietly skip is the defect this project has
    now found three times.
    """

    @pytest.fixture(autouse=True)
    def _needs_the_engine(self):
        pytest.importorskip(
            "pdbfixer", reason="applying a mutation requires the [md] extra")
        pytest.importorskip(
            "openmm", reason="applying a mutation requires the [md] extra")

    def test_the_residue_is_replaced(self, structure, tmp_path):
        out = tmp_path / "out.pdb"
        fix_pdb_with_pdbfixer(str(structure), str(out), mutations=["L2A"])
        assert " ALA A   2" in out.read_text(encoding="utf-8")

    def test_a_mismatch_stops_the_run(self, structure, tmp_path):
        """The guard this feature exists for.

        Applying it anyway would produce a structure, a trajectory and a
        set of results, all describing a protein nobody chose.
        """
        with pytest.raises(ValueError, match="that position holds LEU"):
            fix_pdb_with_pdbfixer(
                str(structure), str(tmp_path / "o.pdb"), mutations=["V2A"])

    def test_a_residue_that_is_not_there(self, structure, tmp_path):
        with pytest.raises(ValueError, match="does not contain"):
            fix_pdb_with_pdbfixer(
                str(structure), str(tmp_path / "o.pdb"), mutations=["L99A"])

    def test_a_chain_that_is_not_there(self, structure, tmp_path):
        with pytest.raises(ValueError, match="chain"):
            fix_pdb_with_pdbfixer(
                str(structure), str(tmp_path / "o.pdb"),
                mutations=["L2A"], mutation_chain="Z")

    def test_no_mutations_changes_nothing_about_the_path(
            self, structure, tmp_path):
        """The feature is inert when unused.

        Every existing study passes no mutations, so the guard and the
        rewrite must not run at all rather than run and decide to do
        nothing.
        """
        out = tmp_path / "out.pdb"
        fix_pdb_with_pdbfixer(str(structure), str(out))
        assert " LEU A   2" in out.read_text(encoding="utf-8")


class TestItReachesAStudyFile:
    def test_the_schema_offers_it(self):
        from fastmdxplora.config.schema import PHASE_SCHEMAS

        names = {f.name for f in PHASE_SCHEMAS["setup"].fields}
        assert {"mutations", "mutation_chain"} <= names

    def test_the_help_says_the_side_chain_is_placed_not_modelled(self):
        """Because it changes what equilibration is doing.

        A deposited side chain sits where crystallography put it; a
        replacement is placed geometrically, and a reader choosing this
        setting is owed that difference.
        """
        from fastmdxplora.config.schema import PHASE_SCHEMAS

        field = PHASE_SCHEMAS["setup"].get("mutations")
        assert "geometrically" in field.help
        assert "checked against the structure" in field.help
