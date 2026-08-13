"""`ligand_pka`: PROPKA's answer, filtered to the copy that was asked about.

The function runs PROPKA on the whole complex and then keeps only the groups
belonging to one residue -- this ligand, this chain, this number -- because a
protein has hundreds of ionizable groups and the question is about one
molecule's. The filtering, the backbone pseudo-group exclusion, the captured
warnings about half-modelled residues, and the two failure modes are all
logic that runs identically whatever PROPKA said, so PROPKA is played here by
a stand-in whose answers are chosen.
"""

from __future__ import annotations

import logging
from pathlib import Path
from types import SimpleNamespace

import pytest

from fastmdxplora.setup.protonation import (
    ProtonationError,
    ligand_pka,
)
import fastmdxplora.setup.protonation as protonation_module


def _atom(resname="BNZ", chain="A", resseq=400):
    return SimpleNamespace(res_name=resname, chain_id=chain, res_num=resseq)


def _group(atom, group_type="COO", pka=4.2, model=3.8):
    return SimpleNamespace(atom=atom, type=group_type,
                           pka_value=pka, model_pka=model)


def _molecule(groups):
    conformation = SimpleNamespace(groups=list(groups))
    return SimpleNamespace(conformation_names=["1A"],
                           conformations={"1A": conformation})


def _propka_saying(molecule=None, *, raises=None, warns=()):
    """A stand-in for propka.run: `single` returns what it is told to."""

    def single(path, optargs=(), write_pka=False):
        for message in warns:
            logging.getLogger("propka").warning(message)
        if raises:
            raise raises
        return molecule

    return SimpleNamespace(single=single)


@pytest.fixture
def fake_propka(monkeypatch):
    def install(**kwargs):
        monkeypatch.setattr(protonation_module, "_propka",
                            lambda: _propka_saying(**kwargs))
    return install


class TestTheFilter:
    def test_only_the_asked_for_copy_is_returned(self, fake_propka) -> None:
        """The complex is full of ionizable groups; the answer is one
        residue's."""
        fake_propka(molecule=_molecule([
            _group(_atom("ASP", "A", 20), pka=3.9),      # the protein's
            _group(_atom("BNZ", "B", 400), pka=5.0),     # the other copy's
            _group(_atom("BNZ", "A", 401), pka=5.1),     # the wrong number
            _group(_atom("BNZ", "A", 400), pka=4.6),     # the one asked for
        ]))
        found = ligand_pka(Path("complex.pdb"), "BNZ", chain="A", resseq=400)
        assert [g.pka for g in found] == [4.6]
        assert found[0].resname == "BNZ" and found[0].chain == "A"

    def test_case_of_the_residue_name_does_not_matter(self, fake_propka) -> None:
        fake_propka(molecule=_molecule([_group(_atom("BNZ", "A", 400))]))
        assert ligand_pka(Path("c.pdb"), "bnz", chain="A", resseq=400)

    def test_a_backbone_pseudo_group_is_not_titratable(self, fake_propka) -> None:
        """PROPKA reports backbone groups with both pKas at zero; treating
        one as ionizable would invent a titration."""
        fake_propka(molecule=_molecule([
            _group(_atom(), group_type="BBN", pka=0.0, model=0.0),
            _group(_atom(), pka=4.6),
        ]))
        found = ligand_pka(Path("c.pdb"), "BNZ", chain="A", resseq=400)
        assert [g.group_type for g in found] == ["COO"]

    def test_no_ionizable_group_is_an_empty_list_not_an_error(
            self, fake_propka) -> None:
        """For a ligand like benzene, nothing is the correct answer."""
        fake_propka(molecule=_molecule([]))
        assert ligand_pka(Path("c.pdb"), "BNZ", chain="A", resseq=400) == []


class TestWhatPropkaComplainsAbout:
    def test_half_modelled_residues_are_counted_and_logged(
            self, fake_propka, caplog) -> None:
        """Surface side chains with poor density are ordinary; the run notes
        how many there were rather than interleaving PROPKA's complaints with
        the phase output."""
        fake_propka(
            molecule=_molecule([_group(_atom())]),
            warns=("Unexpected number of atoms in residue GLU 12 A",
                   "Unexpected number of atoms in residue LYS 90 A",
                   "something else entirely"),
        )
        with caplog.at_level(logging.INFO):
            ligand_pka(Path("c.pdb"), "BNZ", chain="A", resseq=400)
        said = " ".join(record.getMessage() for record in caplog.records)
        assert "2 residues" in said and "missing side-chain atoms" in said

    def test_propkas_own_logging_level_is_restored(self, fake_propka) -> None:
        propka_logger = logging.getLogger("propka")
        before = propka_logger.level
        fake_propka(molecule=_molecule([]))
        ligand_pka(Path("c.pdb"), "BNZ", chain="A", resseq=400)
        assert propka_logger.level == before
        assert not propka_logger.handlers or all(
            not isinstance(h, protonation_module._CapturingHandler)
            for h in propka_logger.handlers)


class TestFailures:
    def test_a_propka_failure_names_the_ligand(self, fake_propka) -> None:
        fake_propka(raises=RuntimeError("unparseable coordinates"))
        with pytest.raises(ProtonationError, match="BNZ"):
            ligand_pka(Path("c.pdb"), "BNZ", chain="A", resseq=400)

    def test_the_level_is_restored_even_on_failure(self, fake_propka) -> None:
        propka_logger = logging.getLogger("propka")
        before = propka_logger.level
        fake_propka(raises=RuntimeError("boom"))
        with pytest.raises(ProtonationError):
            ligand_pka(Path("c.pdb"), "BNZ", chain="A", resseq=400)
        assert propka_logger.level == before

    def test_no_conformation_is_its_own_error(self, fake_propka) -> None:
        fake_propka(molecule=SimpleNamespace(conformation_names=[],
                                             conformations={}))
        with pytest.raises(ProtonationError, match="no conformation"):
            ligand_pka(Path("c.pdb"), "BNZ", chain="A", resseq=400)
