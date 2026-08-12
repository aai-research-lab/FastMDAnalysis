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
