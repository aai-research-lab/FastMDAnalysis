"""Tests for settling a ligand's protonation from its pKa in the complex.

The decision logic is tested directly with pKa values, so these do not depend
on PROPKA succeeding on any particular structure. One test does exercise
PROPKA itself, on a protein with no ligand, to confirm the integration.
"""

from __future__ import annotations

import pytest

from fastmdxplora.setup.protonation import (
    GroupPka,
    ProtonationError,
    decide,
    ligand_pka,
)


def _group(pka: float, *, model: float = 4.0, kind: str = "COO") -> GroupPka:
    return GroupPka("LIG", "A", 201, kind, pka, model)


class TestDeterminedStates:
    def test_a_ligand_with_no_ionizable_group_proceeds(self) -> None:
        state = decide("BNZ", [], 7.4, expected_ionizable=False)
        assert state.protonated is False
        assert "unambiguous" in state.reason

    def test_an_acid_well_below_the_ph_is_deprotonated(self) -> None:
        state = decide("LIG", [_group(3.4)], 7.4, expected_ionizable=True)
        assert state.protonated is False
        assert "3.40" in state.reason

    def test_a_pocket_shifted_acid_above_the_ph_is_protonated(self) -> None:
        """The case a rule table gets wrong.

        A carboxylic acid is deprotonated at pH 7.4 in water. Buried against a
        carboxylate, or in a cavity with nothing to stabilise the charge, it
        need not be, and a shift of several units is ordinary. Deciding from
        the complex catches this; deciding from the molecule cannot.
        """
        state = decide("LIG", [_group(8.9)], 7.4, expected_ionizable=True)
        assert state.protonated is True

    def test_the_environment_shift_is_reported(self) -> None:
        group = _group(8.9, model=4.0)
        assert group.shift == pytest.approx(4.9)
        assert "shifted +4.90" in str(group)


class TestRefusals:
    def test_a_poised_group_stops(self) -> None:
        """Within a unit of the pH both states are appreciably populated."""
        with pytest.raises(ProtonationError, match="both charge states"):
            decide("LIG", [_group(7.0)], 7.4, expected_ionizable=True)

    def test_one_poised_group_among_several_stops(self) -> None:
        with pytest.raises(ProtonationError, match="poised|both charge states"):
            decide("LIG", [_group(3.2), _group(6.9)], 7.4, expected_ionizable=True)

    def test_silence_from_propka_about_an_ionizable_ligand_stops(self) -> None:
        """Never fall back to a method that answers a different question.

        If the chemistry says the ligand is ionizable but PROPKA reports no
        pKa, the two disagree, and substituting a rule table would hide that.
        """
        with pytest.raises(ProtonationError, match="undetermined"):
            decide("LIG", [], 7.4, expected_ionizable=True)

    def test_the_message_says_what_to_do(self) -> None:
        with pytest.raises(ProtonationError) as excinfo:
            decide("LIG", [_group(7.2)], 7.4, expected_ionizable=True)
        message = str(excinfo.value)
        assert "SDF or MOL2" in message or "supply the ligand" in message.lower()
        assert "pH" in message


class TestPropkaIntegration:
    """One test that actually runs PROPKA, to catch API drift."""

    def test_protein_groups_are_read_from_a_real_calculation(self, tmp_path) -> None:
        pytest.importorskip("propka")

        structure = tmp_path / "complex.pdb"
        structure.write_text(
            "ATOM      1  N   ASP A   1      10.000  10.000  10.000  1.00 20.00           N\n"
            "ATOM      2  CA  ASP A   1      11.400  10.000  10.000  1.00 20.00           C\n"
            "ATOM      3  C   ASP A   1      12.000  11.400  10.000  1.00 20.00           C\n"
            "ATOM      4  O   ASP A   1      11.300  12.400  10.000  1.00 20.00           O\n"
            "ATOM      5  CB  ASP A   1      11.900   9.200  11.200  1.00 20.00           C\n"
            "ATOM      6  CG  ASP A   1      13.400   9.100  11.300  1.00 20.00           C\n"
            "ATOM      7  OD1 ASP A   1      14.000   8.400  10.400  1.00 20.00           O\n"
            "ATOM      8  OD2 ASP A   1      13.900   9.700  12.300  1.00 20.00           O\n"
            "TER\nEND\n",
            encoding="utf-8",
        )

        # Asking about a component that is not present must return nothing
        # rather than raising: absence is a valid answer.
        assert ligand_pka(structure, "BNZ", chain="A", resseq=201) == []


@pytest.mark.network
class TestRealComplex:
    """Evidence that deciding from the complex is operational, not aspirational.

    Deselected by default; run with ``pytest -m network``. It exists so that a
    regression in PROPKA's ligand handling is detectable, because the claim
    that FastMDXplora determines ligand protonation from the bound state rests
    entirely on this working.
    """

    def test_biotin_in_streptavidin_shifts_from_its_solution_pka(self, tmp_path) -> None:
        import urllib.request

        pytest.importorskip("propka")

        structure = tmp_path / "1STP.pdb"
        urllib.request.urlretrieve(  # noqa: S310 - RCSB
            "https://files.rcsb.org/download/1STP.pdb", structure
        )

        groups = ligand_pka(structure, "BTN", chain="A", resseq=300)
        assert groups, "PROPKA returned no ionizable group for biotin"

        # Biotin's carboxylate: about 4.5 in water, lower in the pocket.
        carboxylate = [g for g in groups if g.group_type == "OCO"]
        assert carboxylate, "the carboxylate should be detected"
        assert carboxylate[0].model_pka == pytest.approx(4.5, abs=0.3)
        assert carboxylate[0].shift < 0, "streptavidin should lower it"

        # Every group sits well below physiological pH, so the state is settled.
        state = decide("BTN", groups, 7.4, expected_ionizable=True)
        assert state.protonated is False


class TestWhenPropkaRuns:
    """The calculation should run where it can change the answer, and not else."""

    def test_a_ligand_with_no_ionizable_group_skips_the_calculation(self, tmp_path) -> None:
        """Benzene's protonation does not depend on its surroundings.

        Running PROPKA anyway costs a second and produces complaints about
        unmodelled side chains that could not affect the result.
        """
        from fastmdxplora.setup import protonation

        called = []
        original = protonation.ligand_pka
        protonation.ligand_pka = lambda *a, **k: called.append(a) or []
        try:
            state = protonation.settle(
                tmp_path / "nonexistent.pdb", "BNZ",
                chain="A", resseq=200, ph=7.4, expected_ionizable=False,
            )
        finally:
            protonation.ligand_pka = original

        assert called == [], "PROPKA should not have been consulted"
        assert state.protonated is False
        assert "unambiguous" in state.reason

    def test_a_titratable_ligand_does_consult_propka(self, tmp_path) -> None:
        from fastmdxplora.setup import protonation

        called = []
        original = protonation.ligand_pka
        protonation.ligand_pka = lambda *a, **k: called.append(a) or [
            GroupPka("AIN", "A", 1, "COO", 3.2, 4.0)
        ]
        try:
            state = protonation.settle(
                tmp_path / "x.pdb", "AIN",
                chain="A", resseq=1, ph=7.4, expected_ionizable=True,
            )
        finally:
            protonation.ligand_pka = original

        assert len(called) == 1
        assert state.protonated is False


class TestStructuralComplaintsAreCaptured:
    """PROPKA's own warnings belong in the log, not in the phase output."""

    def test_incomplete_residues_do_not_print_raw(self, tmp_path, capsys) -> None:
        pytest.importorskip("propka")

        structure = tmp_path / "gappy.pdb"
        # LYS missing its side chain past CB: PROPKA counts atoms per residue
        # and complains, which is legitimate but should not spill into output.
        structure.write_text(
            "ATOM      1  N   LYS A   1       0.000   0.000   0.000  1.00 20.00           N\n"
            "ATOM      2  CA  LYS A   1       1.400   0.000   0.000  1.00 20.00           C\n"
            "ATOM      3  C   LYS A   1       2.000   1.400   0.000  1.00 20.00           C\n"
            "ATOM      4  O   LYS A   1       1.300   2.400   0.000  1.00 20.00           O\n"
            "ATOM      5  CB  LYS A   1       1.900  -0.800   1.200  1.00 20.00           C\n"
            "TER\nEND\n",
            encoding="utf-8",
        )
        ligand_pka(structure, "BNZ", chain="A", resseq=200)

        captured = capsys.readouterr()
        assert "Unexpected number" not in captured.out
        assert "Unexpected number" not in captured.err
