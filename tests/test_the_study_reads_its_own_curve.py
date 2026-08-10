"""The numbers a reader wants, computed by the study rather than by hand.

A finished umbrella study wrote a curve and left the reading of it to
whoever opened the file. Reading it by hand is how it gets misread: the grid
spans wherever the coordinate went, the windows covered a part of it, and a
minimum taken across the whole grid references a region nothing visited. A
real study looked to have a 164 kJ/mol barrier that way, against 11 measured
where the windows actually were.

And a torsion is a circle. Nine windows from -30 to 165 degrees measure one
of the two paths between their minima; the other is not a boundary, it is
195 degrees nobody put a window on.
"""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

from fastmdxplora.simulation.umbrella import (
    describe_pmf, warn_if_a_circle_is_left_open)


def _curve(lo_deg=-30.0, hi_deg=165.0, height=11.0):
    """A grid over the whole turn, sampled only between lo and hi."""
    c = np.radians(np.linspace(-180, 180, 60))
    inside = (c >= np.radians(lo_deg)) & (c <= np.radians(hi_deg))
    x = np.clip((c - np.radians(lo_deg)) /
                np.radians(hi_deg - lo_deg), 0, 1)
    e = np.where(inside, height * np.sin(x * np.pi) ** 2, np.nan)
    return c, e


class TestTheBarrierIsMeasuredWhereTheWindowsWere:
    def test_it_reports_a_barrier_and_where(self) -> None:
        c, e = _curve()
        got = describe_pmf(c, e, (np.radians(-30), np.radians(165)))
        assert got["barrier_kjmol"] == pytest.approx(11.0, abs=0.5)
        assert -30 < np.degrees(got["barrier_at"]) < 165

    def test_the_covered_range_is_stated(self) -> None:
        """So a consumer cannot reference against absence without being told
        the range exists."""
        c, e = _curve()
        got = describe_pmf(c, e, (np.radians(-30), np.radians(165)))
        lo, hi = np.degrees(got["covered"])
        assert -35 < lo < -25 and 160 < hi < 170

    def test_the_minima_are_found(self) -> None:
        c = np.radians(np.linspace(-30, 165, 40))
        e = 6.0 * np.cos(np.linspace(0, 4 * np.pi, 40)) + 6.0
        got = describe_pmf(c, e)
        assert len(got["minima"]) >= 1
        assert got["minima"][0]["free_energy_kjmol"] <= \
            got["minima"][-1]["free_energy_kjmol"]

    def test_a_curve_of_gaps_reports_nothing(self) -> None:
        """Rather than a barrier of nan, or of zero."""
        c = np.radians(np.linspace(-180, 180, 20))
        got = describe_pmf(c, np.full(20, np.nan))
        assert got["barrier_kjmol"] is None and got["minima"] == []

    def test_an_empty_curve_is_not_an_error(self) -> None:
        assert describe_pmf([], [])["barrier_kjmol"] is None


class TestAPeriodicCoordinateHasNoEnds:
    @staticmethod
    def _plan(centres_deg, variable="torsion"):
        return SimpleNamespace(
            collective_variable=variable,
            windows=[SimpleNamespace(centre=float(np.radians(c)))
                     for c in centres_deg])

    def test_a_partly_tiled_torsion_is_flagged(self) -> None:
        said = warn_if_a_circle_is_left_open(
            self._plan(np.linspace(-30, 165, 9)))
        assert said is not None
        assert "165 degrees with no window" in said

    def test_a_fully_tiled_turn_is_not(self) -> None:
        said = warn_if_a_circle_is_left_open(
            self._plan(np.linspace(-180, 150, 12)))
        assert said is None

    def test_a_distance_is_never_a_circle(self) -> None:
        """A distance has ends. Warning about its 'uncovered arc' would be
        nonsense, and would train people to ignore the message."""
        assert warn_if_a_circle_is_left_open(
            self._plan(np.linspace(-30, 165, 9), "distance")) is None

    def test_one_window_says_nothing(self) -> None:
        assert warn_if_a_circle_is_left_open(self._plan([0.0])) is None

    def test_the_message_says_what_to_do(self) -> None:
        said = warn_if_a_circle_is_left_open(
            self._plan(np.linspace(-30, 165, 9)))
        assert "Tile the full turn" in said
