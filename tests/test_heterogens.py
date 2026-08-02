"""Tests for deciding what a structure's heterogens mean.

The contract under test is narrow and important: a component is simulated
only when the structure determines that it should be, discarded only when the
structure determines that it should be, and otherwise setup refuses. A wrong
answer here produces a trajectory that looks correct and is not.
"""

from __future__ import annotations

import tempfile
from pathlib import Path

import pytest

from fastmdxplora.setup.heterogens import (
    Action,
    AmbiguousStructureError,
    resolve,
)


def _atom(record, serial, name, altloc, resname, chain, seq, x, y, z,
          occupancy=1.00, element="C"):
    return (
        f"{record:<6}{serial:>5} {name:<4}{altloc:1}{resname:>3} "
        f"{chain:1}{seq:>4}    {x:8.3f}{y:8.3f}{z:8.3f}"
        f"{occupancy:6.2f}  0.00          {element:>2}"
    )


def _structure(lines) -> Path:
    path = Path(tempfile.mkdtemp()) / "structure.pdb"
    path.write_text("\n".join(lines) + "\nEND\n", encoding="utf-8")
    return path


PROTEIN = [
    _atom("ATOM", 1, "N", " ", "ALA", "A", 1, 0, 0, 0, element="N"),
    _atom("ATOM", 2, "CA", " ", "ALA", "A", 1, 1.4, 0, 0),
    _atom("ATOM", 3, "O", " ", "ALA", "A", 1, 2.4, 0.5, 0, element="O"),
]


def _actions(path, **kwargs):
    return {d.resname: d.action for d in resolve(path, **kwargs)}


class TestDeterminateStructures:
    """Cases where the structure does answer the question."""

    def test_ligand_kept_buffer_and_water_dropped(self) -> None:
        """The 4W52 case: benzene is the point, HEPES and water are not."""
        lines = list(PROTEIN)
        lines += [_atom("HETATM", 10 + i, f"C{i}", " ", "BNZ", "A", 201,
                        5 + i * 0.1, 5, 5) for i in range(6)]
        lines += [_atom("HETATM", 30 + i, f"C{i}", " ", "EPE", "A", 301,
                        20 + i * 0.1, 20, 20) for i in range(8)]
        lines += [_atom("HETATM", 50, "O", " ", "HOH", "A", 401, 9, 9, 9,
                        element="O")]

        actions = _actions(_structure(lines))
        assert actions["BNZ"] is Action.SIMULATE
        assert actions["EPE"] is Action.DISCARD
        assert actions["HOH"] is Action.DISCARD

    def test_coordinated_metal_kept_free_metal_dropped(self) -> None:
        """A catalytic zinc is structural; a stray sodium is liquor."""
        lines = [
            _atom("ATOM", 1, "NE2", " ", "HIS", "A", 1, 0, 0, 0, element="N"),
            _atom("HETATM", 2, "ZN", " ", "ZN", "A", 101, 2.1, 0, 0, element="ZN"),
            _atom("HETATM", 3, "NA", " ", "NA", "A", 102, 30, 30, 30, element="NA"),
        ]
        actions = _actions(_structure(lines))
        assert actions["ZN"] is Action.SIMULATE
        assert actions["NA"] is Action.DISCARD

    def test_several_distinct_ligands_are_all_simulated(self) -> None:
        lines = list(PROTEIN)
        lines += [_atom("HETATM", 10, "C1", " ", "BNZ", "A", 201, 5, 5, 5)]
        lines += [_atom("HETATM", 11, "C1", " ", "STI", "A", 202, 9, 9, 9)]
        actions = _actions(_structure(lines))
        assert actions["BNZ"] is Action.SIMULATE
        assert actions["STI"] is Action.SIMULATE

    def test_major_conformation_wins_when_occupancy_separates_them(self) -> None:
        lines = list(PROTEIN)
        lines += [
            _atom("HETATM", 10, "C1", "A", "STI", "A", 201, 5, 5, 5, occupancy=0.70),
            _atom("HETATM", 11, "C1", "B", "STI", "A", 201, 6, 5, 5, occupancy=0.30),
        ]
        assert _actions(_structure(lines))["STI"] is Action.SIMULATE

    def test_water_can_be_retained_on_request(self) -> None:
        lines = list(PROTEIN)
        lines += [_atom("HETATM", 50, "O", " ", "HOH", "A", 401, 9, 9, 9,
                        element="O")]
        actions = _actions(_structure(lines), keep_water=True)
        assert actions["HOH"] is Action.SIMULATE

    def test_an_apo_structure_needs_no_decisions(self) -> None:
        assert resolve(_structure(PROTEIN)) == []


class TestRefusals:
    """Cases where guessing would produce a plausible, wrong trajectory."""

    def test_covalent_adduct_stops_setup(self) -> None:
        """A covalent ligand is part of the macromolecule, not a free molecule."""
        lines = list(PROTEIN)
        lines += ["LINK         C1  LIG A 201                 SG  CYS A   1"]
        lines += [_atom("HETATM", 10, "C1", " ", "LIG", "A", 201, 3, 0, 0)]

        with pytest.raises(AmbiguousStructureError, match="covalently bonded"):
            resolve(_structure(lines))

    def test_unknown_component_stops_setup(self) -> None:
        lines = list(PROTEIN)
        lines += [_atom("HETATM", 10, "C1", " ", "UNL", "A", 201, 5, 5, 5)]

        with pytest.raises(AmbiguousStructureError, match="does not say what"):
            resolve(_structure(lines))

    def test_cofactor_needing_special_parameters_stops_setup(self) -> None:
        """Heme is real and essential; OpenFF cannot describe it."""
        lines = list(PROTEIN)
        lines += [_atom("HETATM", 10 + i, f"C{i}", " ", "HEM", "A", 201, 5, 5, 5)
                  for i in range(4)]

        with pytest.raises(AmbiguousStructureError, match="heme"):
            resolve(_structure(lines))

    def test_tied_conformations_stop_setup(self) -> None:
        """Equal occupancy carries no information about which is real."""
        lines = list(PROTEIN)
        lines += [
            _atom("HETATM", 10, "C1", "A", "STI", "A", 201, 5, 5, 5, occupancy=0.50),
            _atom("HETATM", 11, "C1", "B", "STI", "A", 201, 6, 5, 5, occupancy=0.50),
        ]
        with pytest.raises(AmbiguousStructureError, match="occupancy"):
            resolve(_structure(lines))

    def test_partly_coordinated_metal_stops_setup(self) -> None:
        """Keeping or dropping the component as a whole would both be wrong."""
        lines = [
            _atom("ATOM", 1, "NE2", " ", "HIS", "A", 1, 0, 0, 0, element="N"),
            _atom("HETATM", 2, "ZN", " ", "ZN", "A", 101, 2.1, 0, 0, element="ZN"),
            _atom("HETATM", 3, "ZN", " ", "ZN", "A", 102, 40, 40, 40, element="ZN"),
        ]
        with pytest.raises(AmbiguousStructureError, match="coordinated"):
            resolve(_structure(lines))

    def test_the_message_names_everything_unresolved(self) -> None:
        """One run should surface every problem, not the first one only."""
        lines = list(PROTEIN)
        lines += [_atom("HETATM", 10, "C1", " ", "UNL", "A", 201, 5, 5, 5)]
        lines += [_atom("HETATM", 20, "C1", " ", "HEM", "A", 202, 9, 9, 9)]

        with pytest.raises(AmbiguousStructureError) as excinfo:
            resolve(_structure(lines))
        message = str(excinfo.value)
        assert "UNL" in message
        assert "HEM" in message
        assert "Nothing has been simulated" in message


class TestParsing:
    """The classifier is only as good as what it reads out of the file."""

    def test_modified_polymer_residues_are_not_heterogens(self) -> None:
        """Selenomethionine and protonation variants are protein."""
        lines = [
            _atom("ATOM", 1, "N", " ", "MSE", "A", 1, 0, 0, 0, element="N"),
            _atom("ATOM", 2, "N", " ", "HID", "A", 2, 3, 0, 0, element="N"),
            _atom("ATOM", 3, "N", " ", "CYX", "A", 3, 6, 0, 0, element="N"),
        ]
        assert resolve(_structure(lines)) == []

    def test_malformed_lines_do_not_crash_the_run(self) -> None:
        lines = list(PROTEIN) + ["HETATM    x garbage record"]
        resolve(_structure(lines))  # must not raise
