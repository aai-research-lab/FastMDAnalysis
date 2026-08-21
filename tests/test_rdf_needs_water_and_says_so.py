"""`rdf` could not succeed on a default configuration.

Its default pair is the solute against water oxygens; a default run saves
the solute alone. Two defaults asking for different things, so the analysis
failed on every default run -- reproduced on 1L2Y, where the trajectory
loads 311 atoms and the selection finds no water.

The gate for this already existed. `_build_plan` computes `has_water` and
skips anything declaring `requires_water`, and `water_sites` declares it.
`rdf` did not. The whole defect was a missing attribute.

When it *is* asked for explicitly, the refusal now names the setting to
change, which the simulation phase prints hours earlier in a different log.
"""
import mdtraj as md
import numpy as np
import pytest

from fastmdxplora.analysis.rdf import RadialDistribution, _what_would_fix_it


def _solute_only(n_frames=5, n_atoms=12):
    top = md.Topology()
    chain = top.add_chain()
    residue = top.add_residue("ALA", chain)
    for _ in range(n_atoms):
        top.add_atom("CA", md.element.carbon, residue)
    rng = np.random.default_rng(0)
    traj = md.Trajectory(
        (rng.random((n_frames, n_atoms, 3)) * 2.0).astype(np.float32), top)
    traj.unitcell_lengths = np.tile([4.0, 4.0, 4.0], (n_frames, 1))
    traj.unitcell_angles = np.tile([90.0, 90.0, 90.0], (n_frames, 1))
    return traj


def test_rdf_is_gated_on_water_like_water_sites():
    """The attribute the analysis was missing."""
    assert RadialDistribution.requires_water is True


def test_the_refusal_names_the_setting():
    traj = _solute_only()
    with pytest.raises(ValueError) as excinfo:
        RadialDistribution(selection_a="name CA",
                           selection_b="water and name O").compute(traj)
    message = str(excinfo.value)
    assert "matched no atoms" in message
    assert "save_selection: all" in message


@pytest.mark.parametrize("selection", [
    "water and name O", "resname HOH", "solvent", "WATER and name O",
])
def test_solvent_selections_are_recognised(selection):
    assert "save_selection" in _what_would_fix_it(selection, _solute_only())


@pytest.mark.parametrize("selection", [
    "protein and name CA", "resname BEN", "name ZN", "backbone",
])
def test_a_selection_that_is_not_solvent_gets_no_guess(selection):
    """A selection matching nothing is usually a typo. Guessing at a remedy
    would send the reader after the wrong thing."""
    assert _what_would_fix_it(selection, _solute_only()) == ""


def test_no_remedy_offered_when_water_is_present():
    """Then the empty selection has some other cause, and the analysis does
    not know what it is."""
    top = md.Topology()
    chain = top.add_chain()
    water = top.add_residue("HOH", chain)
    top.add_atom("O", md.element.oxygen, water)
    traj = md.Trajectory(np.zeros((2, 1, 3), dtype=np.float32), top)
    traj.unitcell_lengths = np.tile([4.0, 4.0, 4.0], (2, 1))
    traj.unitcell_angles = np.tile([90.0, 90.0, 90.0], (2, 1))
    assert _what_would_fix_it("water and name O", traj) == ""
