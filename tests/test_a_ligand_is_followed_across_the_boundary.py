"""Ligand pose RMSD reported the periodic box instead of the ligand.

`compute` took raw Cartesian displacement, and `superposed` discards the
unit cell, so nothing followed the ligand across a periodic face. Measured
on the 20 ns T4-lysozyme unbound control: 65 frame-to-frame jumps above
2 nm, the largest clustered at 6.58-6.80 nm between frames 10 ps apart, and
a maximum of 9.49 nm in a box smaller than that.

Two earlier attempts imaged the ligand per frame into the receptor's cell.
Both failed, and the failures are why this file is shaped as it is:

* The first rounded fractional coordinates. That is correct only in a
  near-orthogonal cell, and these runs use a rhombic dodecahedron, where it
  picks a longer image than the true minimum for about 30% of random pairs.
  It passed a *cubic* test fixture and moved the real control's maximum
  from 9.49 nm to 9.22 nm, which is no fix at all. Every fixture here is
  triclinic for that reason.

* The second delegated the minimum image to mdtraj, which handles triclinic
  correctly, and still made the jumps worse: 59 to 84. Per-frame imaging
  answers "where is it now" and cannot make a path continuous -- when the
  nearest image flips between frames, the path jumps though the ligand has
  not moved.

Unwrapping is the third approach: each frame takes the image nearest its
predecessor. Consecutive frames are a saving interval apart, so a step of
nearly a whole lattice vector is an image swap and nothing else.

Every test below compares against an independently constructed continuous
trajectory, so a fix that merely moves the numbers cannot pass.
"""
import mdtraj as md
import numpy as np
import pytest

from fastmdxplora.analysis.ligand_rmsd import _followed_across_the_boundary


EDGE = 6.8
#: A rhombic dodecahedron, which is what `box_shape: dodecahedron` builds and
#: what every benchmark run used.
CELL = np.array([[EDGE, 0.0, 0.0],
                 [0.0, EDGE, 0.0],
                 [EDGE / 2, EDGE / 2, EDGE * np.sqrt(2) / 2]])


def _a_wandering_ligand(n_frames=400, step=0.08, seed=1, n_ligand=6):
    """A continuous walk, and the wrapped trajectory a DCD would hold.

    Returns the trajectory as stored, the true unwrapped ligand positions,
    and the two atom selections.
    """
    top = md.Topology()
    chain = top.add_chain()
    protein = top.add_residue("ALA", chain)
    for _ in range(20):
        top.add_atom("CA", md.element.carbon, protein)
    ligand = top.add_residue("BNZ", chain)
    for _ in range(n_ligand):
        top.add_atom("C", md.element.carbon, ligand)

    rng = np.random.default_rng(seed)
    receptor = np.tile(rng.random((20, 3)) + 1.5, (n_frames, 1, 1))
    ring = rng.random((n_ligand, 3)) * 0.15
    walk = np.cumsum(rng.normal(0, step, (n_frames, 1, 3)), axis=0)
    truth = ring[None] + receptor[0].mean(axis=0) + walk

    stored = np.concatenate([receptor, truth], axis=1)
    fractional = stored @ np.linalg.inv(CELL)
    wrapped = (fractional - np.floor(fractional)) @ CELL

    traj = md.Trajectory(wrapped.astype(np.float32), top)
    traj.unitcell_vectors = np.tile(
        CELL, (n_frames, 1, 1)).astype(np.float32)
    return traj, truth, top.select("resname BNZ"), top.select("name CA")


def _rmsd_from_start(xyz):
    displacement = xyz - xyz[0]
    return np.sqrt(np.mean(np.sum(displacement ** 2, axis=2), axis=1))


def test_it_recovers_the_continuous_path():
    """The test the first two attempts could not pass."""
    traj, truth, ligand, receptor = _a_wandering_ligand()
    followed = _followed_across_the_boundary(traj, ligand, receptor)
    assert np.abs(_rmsd_from_start(followed)
                  - _rmsd_from_start(truth)).max() < 1e-5


@pytest.mark.parametrize("step,frames", [(0.05, 400), (0.15, 800), (0.12, 1200)])
def test_it_holds_across_diffusion_rates(step, frames):
    traj, truth, ligand, receptor = _a_wandering_ligand(
        n_frames=frames, step=step, seed=step and int(step * 100))
    followed = _followed_across_the_boundary(traj, ligand, receptor)
    assert np.abs(_rmsd_from_start(followed)
                  - _rmsd_from_start(truth)).max() < 1e-5


def test_the_ligand_does_not_teleport():
    """A jump of a box repeat between adjacent frames is an image swap."""
    traj, _truth, ligand, receptor = _a_wandering_ligand()
    followed = _rmsd_from_start(_followed_across_the_boundary(
        traj, ligand, receptor))
    assert np.abs(np.diff(followed)).max() < 1.0


def test_the_stored_trajectory_really_does_teleport():
    """Without which the assertion above tests nothing."""
    traj, _truth, ligand, _receptor = _a_wandering_ligand()
    raw = _rmsd_from_start(
        np.asarray(traj.xyz[:, ligand, :], dtype=np.float64))
    assert np.abs(np.diff(raw)).max() > 1.0


def test_a_bound_ligand_is_left_where_it_is():
    """The common case must not move. A ligand that never crosses a face
    has one image, and following it is a no-op."""
    traj, truth, ligand, receptor = _a_wandering_ligand(
        n_frames=300, step=0.002, seed=3)
    followed = _rmsd_from_start(_followed_across_the_boundary(
        traj, ligand, receptor))
    assert followed.max() < 0.2
    assert np.abs(followed - _rmsd_from_start(truth)).max() < 1e-6


def test_a_trajectory_with_no_box_is_not_second_guessed():
    """Nothing to follow it across. The raw coordinates are all there is."""
    traj, _truth, ligand, receptor = _a_wandering_ligand(n_frames=50)
    traj.unitcell_vectors = None
    assert _followed_across_the_boundary(traj, ligand, receptor) is None


def test_rounding_fractional_coordinates_would_not_do():
    """The first attempt's method, shown wrong in the cell that matters.

    Kept as a test so the approach is not tried a fourth time.
    """
    rng = np.random.default_rng(0)
    a = rng.random((3000, 3)) @ CELL
    b = rng.random((3000, 3)) @ CELL

    top = md.Topology()
    chain = top.add_chain()
    residue = top.add_residue("X", chain)
    for _ in range(2):
        top.add_atom("C", md.element.carbon, residue)
    traj = md.Trajectory(
        np.stack([a, b], axis=1).astype(np.float32), top)
    traj.unitcell_vectors = np.tile(
        CELL, (3000, 1, 1)).astype(np.float32)

    true_minimum = md.compute_distances(traj, [[0, 1]], periodic=True)[:, 0]
    delta = b - a
    rounded = np.linalg.norm(
        delta - np.round(delta @ np.linalg.inv(CELL)) @ CELL, axis=1)

    assert (rounded > true_minimum + 1e-4).mean() > 0.1
