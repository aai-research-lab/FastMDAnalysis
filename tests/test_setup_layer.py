
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


class TestRemovalReporting:
    """Removing a component and re-adding it are different events."""

    @staticmethod
    def _fix(tmp_path, **kwargs):
        """Run PDBFixer on a tiny structure, capturing what it logged."""
        import io
        import logging

        pytest.importorskip("pdbfixer")
        from fastmdxplora.setup.pdbfix import fix_pdb_with_pdbfixer

        structure = tmp_path / "in.pdb"
        structure.write_text(
            "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N\n"
            "ATOM      2  CA  ALA A   1       1.400   0.000   0.000  1.00  0.00           C\n"
            "ATOM      3  C   ALA A   1       2.000   1.400   0.000  1.00  0.00           C\n"
            "ATOM      4  O   ALA A   1       1.300   2.400   0.000  1.00  0.00           O\n"
            "HETATM   10  C1  BNZ A 200       9.000   9.000   9.000  1.00  0.00           C\n"
            "HETATM   20  C1  EPE A 300      20.000  20.000  20.000  1.00  0.00           C\n"
            "END\n",
            encoding="utf-8",
        )
        captured = io.StringIO()
        handler = logging.StreamHandler(captured)
        handler.setLevel(logging.INFO)
        package_logger = logging.getLogger("fastmdx")
        package_logger.addHandler(handler)
        # A handler alone is not enough: the logger's own level filters
        # records before any handler sees them, so an INFO message about a
        # reinstated component would never arrive.
        previous = package_logger.level
        package_logger.setLevel(logging.INFO)
        try:
            fix_pdb_with_pdbfixer(str(structure), str(tmp_path / "out.pdb"), **kwargs)
        finally:
            package_logger.removeHandler(handler)
            package_logger.setLevel(previous)
        return captured.getvalue()

    def test_dropping_everything_warns_about_everything(self, tmp_path) -> None:
        text = self._fix(tmp_path)
        assert "Removed heterogens" in text
        assert "BNZ" in text
        assert "EPE" in text

    def test_a_judged_component_is_not_second_guessed(self, tmp_path) -> None:
        """The classifier reported EPE as a buffer, with a reason.

        Warning that it "was removed, and might be the ligand you meant to
        simulate" contradicts that, and invites the user to overturn a correct
        decision.
        """
        text = self._fix(tmp_path, reinstated=("BNZ",), explained=("BNZ", "EPE"))
        assert "Removed heterogens" not in text
        assert "re-added with small-molecule parameters" in text

    def test_an_unjudged_component_still_warns(self, tmp_path) -> None:
        """The safety net stays: silence is only earned by having reasoned."""
        text = self._fix(tmp_path, reinstated=("BNZ",), explained=("BNZ",))
        assert "Removed heterogens" in text
        assert "EPE" in text

    def test_a_reinstated_component_is_not_reported_as_lost(self, tmp_path) -> None:
        """Under the auto policy the ligand comes back with parameters.

        Warning that it was removed describes a loss that did not happen, and
        would send a user looking for a ligand that is in the system.
        """
        text = self._fix(tmp_path, reinstated=("BNZ",))

        assert "re-added with small-molecule parameters" in text
        # The buffer really was discarded, so it still warrants the warning.
        assert "Removed heterogens" in text
        assert "EPE" in text
        # ...but benzene must not appear in the removal warning.
        warning_line = next(
            line for line in text.splitlines() if "Removed heterogens" in line
        )
        assert "BNZ" not in warning_line
