"""A periodic profile has to meet itself, and how far it misses is evidence.

The two ends of a full turn are the same place, so their free energies must
be equal. They arrive there from different windows through a chain of joins,
and nothing in the arithmetic forces them to agree -- which makes the
difference a measurement rather than a fault: it is what the study's own
statistics are worth.

Reported, not constrained. Forcing the ends together would make the number
zero and take the information with it. On a well-sampled synthetic profile it
is zero already; a real study of 0.2 ns windows closed to 2 kJ/mol, which is
the honest size of that study's uncertainty and worth printing beside its
barrier.
"""

from __future__ import annotations

import numpy as np
import pytest

from fastmdxplora.simulation.umbrella import closure_gap, describe_pmf


def _turn(bins: int = 60) -> np.ndarray:
    return np.linspace(-np.pi, np.pi, bins)


def _well(x: np.ndarray, height: float = 12.0) -> np.ndarray:
    return height * (1.0 - np.cos(2.0 * x)) / 2.0


class TestItMeasuresTheMiss:
    def test_a_profile_that_closes_reports_zero(self) -> None:
        c = _turn()
        assert closure_gap(c, _well(c)) == pytest.approx(0.0, abs=1e-9)

    @pytest.mark.parametrize("drift", [0.5, 2.0, 179.6])
    def test_a_drift_is_reported_at_its_size(self, drift: float) -> None:
        """179.6 is what the real study gave before the periodic distance
        was fixed: a ramp with no minimum in it."""
        c = _turn()
        energy = _well(c) + np.linspace(0.0, drift, c.size)
        assert closure_gap(c, energy) == pytest.approx(drift, abs=1e-6)

    def test_it_is_a_size_not_a_sign(self) -> None:
        c = _turn()
        up = _well(c) + np.linspace(0.0, 3.0, c.size)
        down = _well(c) + np.linspace(3.0, 0.0, c.size)
        assert closure_gap(c, up) == pytest.approx(closure_gap(c, down))


class TestItOnlyAppliesToACircle:
    def test_a_partial_turn_has_no_closure_to_measure(self) -> None:
        """Two ends a long way short of meeting are not failing to close;
        they are simply not a circle, and a number here would invite somebody
        to read a drift into an ordinary profile."""
        c = np.radians(np.linspace(-30, 165, 40))
        assert closure_gap(c, 11.0 * np.sin(np.linspace(0, np.pi, 40)) ** 2) is None

    def test_a_curve_of_gaps_reports_nothing(self) -> None:
        assert closure_gap(_turn(), np.full(60, np.nan)) is None

    def test_one_finite_bin_is_not_a_profile(self) -> None:
        energy = np.full(60, np.nan)
        energy[10] = 1.0
        assert closure_gap(_turn(), energy) is None


class TestTheStudyReportsIt:
    def test_it_travels_with_the_barrier(self) -> None:
        c = _turn()
        got = describe_pmf(c, _well(c) + np.linspace(0.0, 2.0, c.size),
                           periodic=True)
        assert got["closure_gap_kjmol"] == pytest.approx(2.0, abs=1e-6)
        assert got["barrier_kjmol"] is not None

    def test_a_partial_study_carries_none_rather_than_zero(self) -> None:
        """Zero would read as 'closed perfectly', which is the opposite of
        'this was never a closed turn'."""
        c = np.radians(np.linspace(-30, 165, 40))
        got = describe_pmf(c, 11.0 * np.sin(np.linspace(0, np.pi, 40)) ** 2,
                           periodic=True)
        assert got["closure_gap_kjmol"] is None
