"""Insertion codes are read from whichever file gave the topology.

They were read only from an external ``.pdb`` topology. A PDB loaded as the
trajectory itself, or an mmCIF topology, lost them, and trypsin's GLY 184A
and TYR 184 came out as two residues numbered 184 in one chain.
"""

from __future__ import annotations

import collections
import gzip
from pathlib import Path

import numpy as np
import pytest

md = pytest.importorskip("mdtraj")

DEPOSITED = Path(__file__).parent / "data" / "assemblies"


@pytest.fixture(autouse=True)
def _no_codes_from_earlier_tests(monkeypatch):
    from fastmdxplora.analysis import residues

    monkeypatch.setattr(residues, "_INSERTION_CODES", {})


def _deposited(tmp_path, name):
    path = tmp_path / name.replace(".gz", "")
    path.write_bytes(gzip.decompress((DEPOSITED / name).read_bytes()))
    return path


def _rmsf_labels(traj):
    from fastmdxplora.analysis.rmsf import RMSF

    table = RMSF(start=0).compute(md.join([traj] * 3))
    assert hasattr(table, "columns"), "a structure with insertion codes is a table"
    return list(zip(table["residue"], table["insertion"]))


def test_a_pdb_loaded_as_the_trajectory_keeps_them(tmp_path):
    from fastmdxplora.analysis.loading import load_trajectory

    traj = load_trajectory(str(_deposited(tmp_path, "3PTB.pdb.gz")))
    labels = _rmsf_labels(traj)

    assert not [k for k, n in collections.Counter(labels).items() if n > 1]
    assert {(184, "A"), (184, ""), (188, "A"), (221, "A")} <= set(labels)


def test_an_mmcif_topology_keeps_them(tmp_path):
    pytest.importorskip("openmm")  # MDTraj reads mmCIF through OpenMM
    from fastmdxplora.analysis.loading import load_trajectory

    cif = _deposited(tmp_path, "3PTB.cif.gz")
    structure = md.load(str(cif))
    structure.save_dcd(str(tmp_path / "one.dcd"))
    traj = load_trajectory(str(tmp_path / "one.dcd"), top=str(cif))
    labels = _rmsf_labels(traj)

    assert not [k for k, n in collections.Counter(labels).items() if n > 1]
    assert {(184, "A"), (184, ""), (188, "A"), (221, "A")} <= set(labels)


def test_the_mmcif_reader_takes_the_columns_by_name(tmp_path):
    """Read as the format writes them, quoted values and all, first model
    only, and keyed as MDTraj names the atoms."""
    from fastmdxplora.analysis import residues

    path = tmp_path / "x.cif"
    path.write_text(
        "data_X\n#\nloop_\n"
        "_atom_site.group_PDB\n_atom_site.id\n_atom_site.label_atom_id\n"
        "_atom_site.label_comp_id\n_atom_site.pdbx_PDB_ins_code\n"
        "_atom_site.auth_seq_id\n_atom_site.auth_comp_id\n_atom_site.auth_asym_id\n"
        "_atom_site.pdbx_PDB_model_num\n"
        "ATOM 1 \"O5'\" GLY A 184 GLY A 1\n"
        "ATOM 2 CA TYR ? 184 TYR A 1\n"
        "ATOM 3 CA GLY A 184 GLY A 2\n"
        "#\n", encoding="utf-8")

    assert residues.remember_insertion_codes(path) == 1
    assert residues._INSERTION_CODES == {(1, "GLY", 184): "A"}


def test_a_file_without_codes_records_none(tmp_path):
    from fastmdxplora.analysis import residues

    assert residues.remember_insertion_codes(_deposited(tmp_path, "1AKE.pdb.gz")) == 0
    assert residues._INSERTION_CODES == {}


def test_the_numbers_of_an_inserted_residue_stay_numbers(tmp_path):
    from fastmdxplora.analysis.loading import load_trajectory

    traj = load_trajectory(str(_deposited(tmp_path, "3PTB.pdb.gz")))
    assert all(isinstance(n, (int, np.integer)) for n, _ in _rmsf_labels(traj))


def _glycines(path, numbers):
    """Four glycines of N, CA, C and O, numbered as given: "184A" has a code."""
    lines, serial = [], 1
    for i, given in enumerate(numbers):
        resseq, code = int(given[:3]), given[3:]
        for k, name in enumerate(("N", "CA", "C", "O")):
            lines.append(f"ATOM  {serial:5d}  {name:<3s} GLY A{resseq:4d}{code:1s}   "
                         f"{3.6 * i + 0.4 * k:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00  0.00           {name[0]}\n")
            serial += 1
    path.write_text("".join(lines) + "END\n")
    return path


def test_a_second_structure_read_in_one_process_has_only_its_own_codes(tmp_path):
    """A structure read after one with an insertion code had that code
    wherever an atom had the same serial, name and number: one with none
    read as having insertions, and its 184 was "184A"."""
    from fastmdxplora.analysis.loading import load_trajectory
    from fastmdxplora.analysis.residues import has_insertions, label

    first = _glycines(tmp_path / "a.pdb", ["183", "184A", "184", "185"])
    second = _glycines(tmp_path / "b.pdb", ["183", "184", "185", "186"])

    a = load_trajectory(str(first))
    assert has_insertions(a.topology)
    assert [label(r, qualified=False) for r in a.topology.residues] == [183, "184A", 184, 185]

    b = load_trajectory(str(second))
    assert not has_insertions(b.topology)
    assert [label(r, qualified=False) for r in b.topology.residues] == [183, 184, 185, 186]

    again = load_trajectory(str(first))
    assert [label(r, qualified=False) for r in again.topology.residues] == [183, "184A", 184, 185]
