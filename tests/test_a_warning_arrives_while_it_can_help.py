"""What a run would tell you, said while you can still act on it.

The setup phase says a great deal: which heterogens it discarded, that this
force field wants hard truncation, that a metal in a site will not be held
there by charge alone. All of it is true and all of it arrives once the run
has begun -- by which point the person who would have changed something has
stopped watching, and on a cluster has gone home.

Most of it is decidable from a structure and a set of settings. Said while
somebody is still choosing, the same sentence changes what they do instead of
explaining what already happened.

None of this refuses anything. Every case here runs perfectly well and
produces a result worth doubting, which is exactly why a validator is the
wrong shape for it.
"""

from __future__ import annotations

import pytest

from fastmdxplora.advisories import Advisory, advise


def _settings(**kwargs):
    return kwargs


class TestAMetalInASite:
    def test_a_zinc_is_worth_saying(self) -> None:
        said = advise({"ion_resnames": ["ZN"]}, _settings())
        assert [a.setting for a in said] == ["forcefield"]

    def test_salt_is_not(self) -> None:
        """Sodium and chloride are meant to move. Warning about them would
        fire on every run and teach people to ignore the rest."""
        assert not advise({"ion_resnames": ["NA", "CL"]}, _settings())

    def test_it_says_what_to_do(self) -> None:
        said = advise({"ion_resnames": ["ZN", "CA"]}, _settings())[0]
        assert "restrain" in said.remedy
        assert "ZN" in said.summary and "CA" in said.summary

    def test_it_reads_the_same_list_the_run_does(self) -> None:
        """Two ideas of what counts as a structural metal would drift, and a
        warning at config time that the run does not repeat is worse than
        none."""
        from pathlib import Path
        import fastmdxplora.advisories as advisories

        source = Path(advisories.__file__).read_text(encoding="utf-8")
        assert "from fastmdxplora.setup.prepare import STRUCTURAL_METALS" in source


class TestABoxTooSmallForTheCutoff:
    """Calibrated against a real refusal. A capped alanine at 1.0 and 1.2 nm
    of padding could not carry a 1.0 nm cutoff; 1.6 was the first that
    cleared."""

    @pytest.mark.parametrize("padding,expected", [
        (1.0, True), (1.2, True), (1.6, False), (1.8, False),
    ])
    def test_it_matches_what_actually_happened(self, padding, expected) -> None:
        said = advise({"extents_angstrom": [9, 7, 7]},
                      _settings(solvent_padding_nm=padding,
                                nonbonded_cutoff_nm=1.0))
        assert bool(said) is expected

    def test_a_real_protein_is_quiet(self) -> None:
        assert not advise({"extents_angstrom": [40, 30, 30]},
                          _settings(solvent_padding_nm=1.0))

    def test_it_says_nothing_without_a_padding(self) -> None:
        """An empty form is where somebody is least able to act on advice."""
        assert not advise({"extents_angstrom": [9, 7, 7]}, _settings())


class TestASwitchTheForceFieldDoesNotWant:
    def test_amber_switched_is_worth_saying(self) -> None:
        said = advise({}, _settings(forcefield="amber14",
                                    use_switching_function=True))
        assert [a.setting for a in said] == ["use_switching_function"]

    def test_charmm_switched_is_not(self) -> None:
        """Switching is what CHARMM wants. Warning there would be wrong."""
        assert not advise({}, _settings(forcefield="charmm36",
                                        use_switching_function=True))

    def test_amber_unswitched_is_not(self) -> None:
        assert not advise({}, _settings(forcefield="amber14",
                                        use_switching_function=False))


class TestALigandWithNoChemistry:
    def test_a_local_structure_with_a_ligand(self) -> None:
        said = advise({"ligand_resnames": ["BNZ"], "path": "181L.pdb"},
                      _settings())
        assert [a.setting for a in said] == ["ligand"]

    def test_supplying_one_settles_it(self) -> None:
        assert not advise({"ligand_resnames": ["BNZ"], "path": "181L.pdb"},
                          _settings(ligand="bnz.sdf"))

    def test_it_says_the_coordinates_do_not_matter(self) -> None:
        """The commonest confusion: an ideal SDF's pose is arbitrary, and the
        ligand is placed where the structure has it."""
        said = advise({"ligand_resnames": ["BNZ"], "path": "181L.pdb"},
                      _settings())[0]
        assert "coordinates do not matter" in said.remedy


class TestADensityNeverEquilibrated:
    def test_no_npt_is_worth_saying(self) -> None:
        said = advise({}, _settings(npt_steps=0))
        assert [a.setting for a in said] == ["npt_steps"]

    def test_npt_settles_it(self) -> None:
        assert not advise({}, _settings(npt_steps=50000))


class TestItIsShapedToBeShownBesideTheSetting:
    def test_each_names_the_field_it_is_about(self) -> None:
        """So an interface can put it beside the control rather than in a
        list somewhere else, which is the whole point of saying it early."""
        for said in advise({"ion_resnames": ["ZN"]}, _settings(npt_steps=0)):
            assert said.setting
            assert said.summary and said.detail and said.remedy

    def test_nothing_at_all_is_not_an_error(self) -> None:
        assert advise(None, None) == []

    def test_the_gui_serves_them(self) -> None:
        from pathlib import Path
        import fastmdxplora.gui.server as server

        source = Path(server.__file__).read_text(encoding="utf-8")
        assert 'info["advisories"]' in source

    def test_they_can_be_read_as_one_sentence(self) -> None:
        said = Advisory("x", "A.", "B.", "C.")
        assert said.as_text() == "A. B. C."


class TestARunSaysThemToo:
    """The GUI is not where a cluster user is. They run `fastmdx explore` and
    walk away, so advice they must remember to ask for is advice they will
    not see -- and a separate `check` subcommand is exactly that. The run
    says them itself, before the expensive part.
    """

    @staticmethod
    def _structure(tmp_path, entries):
        lines = []
        for index, (atom, residue) in enumerate(entries, 1):
            lines.append(
                f"HETATM{index:5d}  {atom:<3s} {residue:>3s} A{900 + index:4d}    "
                f"{index * 2:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00  0.00")
        path = tmp_path / "system.pdb"
        path.write_text("\n".join(lines) + "\nEND\n", encoding="utf-8")
        return path

    def _run(self, tmp_path, entries, options, caplog):
        import logging

        from fastmdxplora.orchestrator import FastMDXplora

        app = FastMDXplora.__new__(FastMDXplora)
        app.system = str(self._structure(tmp_path, entries))
        with caplog.at_level(logging.WARNING):
            app._say_what_is_worth_knowing(options)
        return caplog.text

    def test_a_metal_is_named_before_the_run(self, tmp_path, caplog) -> None:
        said = self._run(tmp_path, [("NE2", "HIS"), ("ZN", "ZN")], {}, caplog)
        assert "ZN" in said and "drift out of its site" in said

    def test_salt_alone_says_nothing(self, tmp_path, caplog) -> None:
        said = self._run(tmp_path, [("NE2", "HIS"), ("NA", "NA")], {}, caplog)
        assert "drift out of its site" not in said

    def test_a_setting_is_read_from_the_options(self, tmp_path, caplog) -> None:
        said = self._run(tmp_path, [("NE2", "HIS")],
                         {"simulation": {"npt_steps": 0}}, caplog)
        assert "density solvation gave it" in said

    def test_it_runs_before_the_first_phase(self) -> None:
        """After setup would be after the point of saying it."""
        from pathlib import Path
        import fastmdxplora.orchestrator as orchestrator

        source = Path(orchestrator.__file__).read_text(encoding="utf-8")
        assert (source.index("_say_what_is_worth_knowing(merged_options)")
                < source.index("for phase in plan:"))

    def test_a_structure_that_will_not_read_is_not_an_error(
            self, tmp_path, caplog) -> None:
        import logging

        from fastmdxplora.orchestrator import FastMDXplora

        broken = tmp_path / "system.pdb"
        broken.write_text("not a structure\n", encoding="utf-8")
        app = FastMDXplora.__new__(FastMDXplora)
        app.system = str(broken)
        with caplog.at_level(logging.WARNING):
            app._say_what_is_worth_knowing({})

    def test_no_structure_at_all_is_not_an_error(self, caplog) -> None:
        from fastmdxplora.orchestrator import FastMDXplora

        app = FastMDXplora.__new__(FastMDXplora)
        app.system = None
        app._say_what_is_worth_knowing({"simulation": {"npt_steps": 0}})


class TestTheStructureReportNamesItsIons:
    def test_which_ions_not_only_how_many(self, tmp_path) -> None:
        """A count cannot tell a structural zinc from added salt, and that
        distinction is what decides whether the metal advisory applies."""
        from fastmdxplora.gui.structure_info import count_structure

        path = tmp_path / "s.pdb"
        path.write_text(
            "HETATM    1  ZN   ZN A 901       1.000   0.000   0.000  1.00  0.00\n"
            "HETATM    2  NA   NA A 902       5.000   0.000   0.000  1.00  0.00\n"
            "END\n", encoding="utf-8")
        got = count_structure(path)
        assert got["ion_resnames"] == ["NA", "ZN"]
        assert got["ions"] == 2
