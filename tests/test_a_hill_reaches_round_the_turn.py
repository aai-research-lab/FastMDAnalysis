"""A metadynamics surface on a circle has no wall at the wrap.

The bias is a sum of Gaussians, and the separation was computed straight:
`grid - centre`. On a torsion that is wrong at the boundary. A hill deposited
at +170 degrees is 15 degrees from a grid point at -175, not 345, and computed
straight it contributes nothing there -- so the surface grows an artificial
wall exactly where the coordinate is continuous, and a frame sampled the other
side of it is given a weight that ignores the bias it actually felt.

The same defect as the umbrella stitching, in three more places: the free
energy surface, the bias each frame felt, and the c(t) offset. Found by
looking for it before running a torsion metadynamics study, rather than after.
"""

from __future__ import annotations

import numpy as np
import pytest

from fastmdxplora.analysis.reweighted_averages import felt_bias
from fastmdxplora.simulation.metad_surface import surface_from_hills
from fastmdxplora.simulation.umbrella import PERIODIC_VARIABLES


class _Hills:
    """One hill near the top of the turn."""

    def __init__(self, centre_deg: float = 170.0, sigma: float = 0.35):
        self.centre = np.radians(np.array([centre_deg]))
        self.sigma = np.array([sigma])
        self.height = np.array([1.2])
        self.time_ps = np.array([1.0])
        self.bias_factor = 10.0

    def __len__(self) -> int:
        return 1


class TestTheSurfaceReachesAcrossTheWrap:
    def test_a_point_just_past_the_boundary_feels_the_hill(self) -> None:
        grid = np.radians(np.array([175.0, -175.0]))
        straight = surface_from_hills(_Hills(), grid)
        periodic = surface_from_hills(_Hills(), grid, periodic=True)
        # Computed straight, the two sides differ enormously for points that
        # are ten degrees apart.
        assert abs(straight[0] - straight[1]) > abs(periodic[0] - periodic[1])

    def test_the_surface_is_continuous_around_the_turn(self) -> None:
        """The ends of a full turn are the same place, so a surface built on
        it should not step between them."""
        grid = np.linspace(-np.pi, np.pi, 180)
        got = surface_from_hills(_Hills(centre_deg=170.0), grid, periodic=True)
        assert abs(got[0] - got[-1]) < 0.05

    def test_a_straight_coordinate_is_unchanged(self) -> None:
        """A distance has ends; wrapping one would invent a neighbour."""
        grid = np.linspace(0.5, 2.5, 40)
        hills = _Hills(centre_deg=0.0)
        hills.centre = np.array([1.5])
        assert np.allclose(surface_from_hills(hills, grid),
                           surface_from_hills(hills, grid, periodic=False))


class TestTheFeltBiasReachesToo:
    def _felt(self, value_deg: float, periodic: bool) -> float:
        hills = _Hills()
        return float(felt_bias(
            hills, hills.height, times=np.array([5.0]),
            values=np.radians(np.array([value_deg])), periodic=periodic)[0])

    def test_a_frame_past_the_wrap_feels_the_bias(self) -> None:
        """15 degrees from the hill, and it felt none of it."""
        assert self._felt(-175.0, periodic=False) == pytest.approx(0.0, abs=1e-9)
        assert self._felt(-175.0, periodic=True) > 0.5

    def test_a_frame_beside_the_hill_is_unaffected(self) -> None:
        near = self._felt(175.0, periodic=False)
        assert self._felt(175.0, periodic=True) == pytest.approx(near)


class TestThePeriodicVariablesAreOneList:
    def test_the_wrapping_stays_out_of_the_weighting_maths(self) -> None:
        """`reweight.py` holds generic weighting arithmetic; a periodic
        Gaussian is the sum over images, so the convention is applied by
        summing the plain calculation rather than by teaching it circles."""
        from pathlib import Path
        import fastmdxplora.analysis.reweight as generic

        assert "periodic" not in Path(generic.__file__).read_text(
            encoding="utf-8")

    def test_the_reweighting_uses_the_umbrella_list(self) -> None:
        """Two lists of which coordinates are circular would eventually
        disagree, and the disagreement would be silent."""
        from pathlib import Path
        import fastmdxplora.analysis.reweighted_averages as module

        source = Path(module.__file__).read_text(encoding="utf-8")
        assert "from fastmdxplora.simulation.umbrella import PERIODIC_VARIABLES" \
            in source

    def test_a_torsion_is_a_circle_and_a_bond_angle_is_not(self) -> None:
        """PLUMED's TORSION wraps; its ANGLE has domain [0, pi].

        This asserted that both were circles, which is what the code said
        and not what PLUMED does. The consequence was a surface grid on
        [-pi, pi] for a bond angle -- half of it ground no geometry can
        occupy -- a barrier taken as a maximum across that half, and
        umbrella advice to tile windows onto it. `displacement()` was
        unaffected, being a no-op for values already in [0, pi], so the bias
        and the histogramming were right throughout.
        """
        assert "torsion" in PERIODIC_VARIABLES
        assert "angle" not in PERIODIC_VARIABLES

    def test_a_distance_is_not(self) -> None:
        assert not ({"distance", "radius_of_gyration", "ligand_rmsd"}
                    & PERIODIC_VARIABLES)
