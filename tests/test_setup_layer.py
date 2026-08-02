
import pytest

def test_removed_heterogens_are_reported() -> None:
    """A stripped ligand must not pass silently.

    Removing crystallization additives is usually correct, but a bound ligand
    is indistinguishable from a buffer molecule at this stage. Running 4W52
    without --setup-keep-heterogens discards the benzene and produces what
    looks like a holo run but is apo.
    """
    from fastmdxplora.setup.pdbfix import _heterogen_residue_counts

    class _Residue:
        def __init__(self, name: str) -> None:
            self.name = name

    class _Topology:
        def __init__(self, names: list[str]) -> None:
            self._residues = [_Residue(n) for n in names]

        def residues(self):
            return iter(self._residues)

    # The real composition of PDB 4W52.
    counts = _heterogen_residue_counts(
        _Topology(["ALA", "GLY", "BNZ", "EPE", "HOH", "HOH", "CL", "MSE"])
    )
    assert counts == {"BNZ": 1, "EPE": 1}

    # Copies are counted so the message can say how many were lost.
    assert _heterogen_residue_counts(_Topology(["BNZ", "BNZ"])) == {"BNZ": 2}

    # Nothing to report for a plain protein or nucleic acid.
    assert _heterogen_residue_counts(_Topology(["ALA", "GLY", "HOH"])) == {}
    assert _heterogen_residue_counts(_Topology(["DA", "DC", "DG", "DT"])) == {}


class TestHeterogenRemovalIsAnnounced:
    """Removing a heterogen must never be silent.

    A bound ligand and a buffer molecule look identical in a PDB file at the
    preparation stage. Dropping both without a word produced runs that read as
    holo simulations but were apo.
    """

    @staticmethod
    def _topology(names):
        class Residue:
            def __init__(self, name):
                self.name = name

        class Topology:
            def __init__(self, residue_names):
                self._residues = [Residue(name) for name in residue_names]

            def residues(self):
                return iter(self._residues)

        return Topology(names)

    def test_ligand_and_buffer_are_both_reported(self) -> None:
        from fastmdxplora.setup.pdbfix import _heterogen_residue_counts

        # 4W52: benzene ligand, HEPES buffer, crystallographic waters.
        counts = _heterogen_residue_counts(
            self._topology(["ALA", "LEU"] * 20 + ["BNZ", "EPE"] + ["HOH"] * 80)
        )
        assert counts == {"BNZ": 1, "EPE": 1}

    def test_apo_structure_reports_nothing(self) -> None:
        """Water and common ions are not worth warning about."""
        from fastmdxplora.setup.pdbfix import _heterogen_residue_counts

        counts = _heterogen_residue_counts(
            self._topology(["ALA"] * 30 + ["HOH"] * 100 + ["NA", "CL", "MG"])
        )
        assert counts == {}

    def test_multiple_copies_are_counted(self) -> None:
        from fastmdxplora.setup.pdbfix import _heterogen_residue_counts

        counts = _heterogen_residue_counts(self._topology(["ALA"] * 10 + ["BNZ"] * 3))
        assert counts == {"BNZ": 3}

    def test_modified_residues_are_not_flagged(self) -> None:
        """MSE and protonation variants are protein, not heterogens."""
        from fastmdxplora.setup.pdbfix import _heterogen_residue_counts

        counts = _heterogen_residue_counts(
            self._topology(["ALA", "MSE", "HID", "HIE", "CYX", "DA", "DT"])
        )
        assert counts == {}


class TestMultipleLigands:
    """Several ligands, cofactors, or copies can be parameterized together.

    They must stay distinguishable: two ligands sharing a residue name cannot
    be told apart in the topology, and every analysis that selects by residue
    would silently conflate them.
    """

    def test_one_ligand_uses_the_name_given(self) -> None:
        from fastmdxplora.setup.prepare import _resolve_ligand_names

        assert _resolve_ligand_names(["bnz.sdf"], "BNZ") == ["BNZ"]
        assert _resolve_ligand_names(["x.sdf"], None) == ["LIG"]

    def test_several_ligands_take_their_names_from_the_files(self) -> None:
        from fastmdxplora.setup.prepare import _resolve_ligand_names

        names = _resolve_ligand_names(["BNZ.sdf", "ATP.sdf", "GOL.sdf"], None)
        assert names == ["BNZ", "ATP", "GOL"]

    def test_explicit_names_are_accepted_one_per_ligand(self) -> None:
        from fastmdxplora.setup.prepare import _resolve_ligand_names

        assert _resolve_ligand_names(["a", "b"], ["L01", "L02"]) == ["L01", "L02"]

    def test_one_name_for_several_ligands_is_refused(self) -> None:
        from fastmdxplora.setup.prepare import _resolve_ligand_names

        with pytest.raises(ValueError, match="single ligand name"):
            _resolve_ligand_names(["a", "b"], "LIG")

    def test_duplicate_names_are_refused(self) -> None:
        from fastmdxplora.setup.prepare import _resolve_ligand_names

        with pytest.raises(ValueError, match="must be distinct"):
            _resolve_ligand_names(["a", "b"], ["X", "X"])

    def test_a_mismatched_number_of_names_is_refused(self) -> None:
        from fastmdxplora.setup.prepare import _resolve_ligand_names

        with pytest.raises(ValueError, match="ligand names given"):
            _resolve_ligand_names(["a", "b"], ["X"])

    def test_charges_are_broadcast_or_given_per_ligand(self) -> None:
        from fastmdxplora.setup.prepare import _resolve_ligand_charges

        assert _resolve_ligand_charges(["a", "b"], None) == [None, None]
        assert _resolve_ligand_charges(["a", "b"], [0, -2]) == [0, -2]

        with pytest.raises(ValueError, match="ligand charges given"):
            _resolve_ligand_charges(["a", "b"], [0])
