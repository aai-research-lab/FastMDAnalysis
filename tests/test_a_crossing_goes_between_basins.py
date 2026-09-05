"""A recrossing is travel between the two states, not past a line.

The thresholds were placed a quarter of the way in from the extremes of the
hill range -- a statement about the grid, not about the system. On a full turn
that put them at -90 and +90 degrees, and a real study whose minima sat at -17
and 155 had one of them *inside* the dead band. What was counted was
excursions past an arbitrary line.

It reported 925 against a required 4, so that run's conclusion held. On a
system whose minima sit differently it would not have, and it would not have
said so.

The fourth place a straight-line assumption was made about a circular
coordinate, after the umbrella stitching, the metadynamics surface and the
felt bias.
"""

from __future__ import annotations

import numpy as np

from fastmdxplora.simulation.metad_surface import basins, transitions


def _double_well(bins: int = 61) -> tuple[np.ndarray, np.ndarray]:
    """Minima at 0 and at the wrap."""
    grid = np.radians(np.linspace(-180, 180, bins))
    return grid, 12.0 * (1.0 - np.cos(2.0 * grid)) / 2.0


class TestTheStatesAreFound:
    def test_both_minima_including_the_one_at_the_wrap(self) -> None:
        grid, surface = _double_well()
        found = basins(grid, surface, periodic=True)
        assert found is not None
        degrees = sorted(round(float(np.degrees(x))) for x in found)
        assert degrees == [-180, 0]

    def test_a_duplicated_endpoint_does_not_hide_one(self) -> None:
        """A grid spanning -pi to pi inclusive names the same place twice, so
        the wrap-around neighbour of the first point is the second to last."""
        grid, surface = _double_well()
        assert abs((grid[-1] - grid[0]) - 2 * np.pi) < 1e-9
        assert basins(grid, surface, periodic=True) is not None

    def test_one_minimum_is_no_basins(self) -> None:
        """A coordinate with a single well has nothing to cross to, and a
        count of crossings would be a count of nothing."""
        grid = np.radians(np.linspace(-180, 180, 61))
        assert basins(grid, np.abs(grid), periodic=True) is None

    def test_the_deepest_two_are_chosen(self) -> None:
        grid = np.radians(np.linspace(-180, 180, 181))
        surface = (10.0 * (1 - np.cos(2 * grid)) / 2
                   + 3.0 * np.sin(grid))          # breaks the degeneracy
        found = basins(grid, surface, periodic=True)
        assert found is not None and len(found) == 2


class TestTravelIsCounted:
    def _walk(self, *stops_deg: float) -> np.ndarray:
        return np.radians(np.array([s for stop in stops_deg
                                    for s in [stop] * 5], dtype=float))

    def test_each_journey_counts_once(self) -> None:
        grid, surface = _double_well()
        two = basins(grid, surface, periodic=True)
        walk = self._walk(0, 180, 0, 180)
        assert transitions(walk, between=two, periodic=True) == 3

    def test_the_two_ends_of_the_turn_are_one_basin(self) -> None:
        """-180 and +180 are the same place. Counting a step between them as
        a transition would invent journeys that never happened."""
        grid, surface = _double_well()
        two = basins(grid, surface, periodic=True)
        assert transitions(self._walk(180, -180, 180), between=two,
                           periodic=True) == 0

    def test_rattling_on_the_barrier_counts_nothing(self) -> None:
        grid, surface = _double_well()
        two = basins(grid, surface, periodic=True)
        rattle = np.radians(np.array([88.0, 92.0] * 20))
        assert transitions(rattle, between=two, periodic=True) == 0

    def test_staying_put_counts_nothing(self) -> None:
        grid, surface = _double_well()
        two = basins(grid, surface, periodic=True)
        assert transitions(self._walk(0, 0, 0), between=two,
                           periodic=True) == 0

    def test_an_empty_trajectory_is_zero(self) -> None:
        assert transitions(np.array([]), between=(0.0, np.pi)) == 0

    def test_two_basins_at_the_same_place_are_not_a_journey(self) -> None:
        assert transitions(np.radians(np.array([0.0, 1.0])),
                           between=(0.5, 0.5)) == 0


class TestTheRecordSaysWhichMeasureItUsed:
    def test_the_definition_names_the_basins(self) -> None:
        """The two measures answer different questions, and the number alone
        does not distinguish them."""
        from pathlib import Path
        import fastmdxplora.simulation.metad_surface as module

        source = Path(module.__file__).read_text(encoding="utf-8")
        assert "travel between the basins at" in source
        assert "there is no second basin to travel to" in source

    def test_where_the_basins_were_is_recorded(self) -> None:
        from pathlib import Path
        import fastmdxplora.simulation.metad_surface as module

        source = Path(module.__file__).read_text(encoding="utf-8")
        assert '"basins":' in source
