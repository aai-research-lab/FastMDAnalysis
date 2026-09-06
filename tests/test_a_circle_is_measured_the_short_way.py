"""A torsion is a circle, and the restraint has to know it.

The bias on a window was `0.5 * k * (x - centre)**2`, a straight-line
subtraction. On a circle that is wrong at the wrap: a sample at +170 degrees
is ten degrees from a window held at -180, not three hundred and fifty. The
code charged it 933 kJ/mol instead of 0.76, which is a Boltzmann weight wrong
by a factor of 10^162.

A real study found it. Twelve windows tiling a full turn of alanine
dipeptide's psi returned a monotonic ramp of 180 kJ/mol with no minimum
anywhere in it -- while every check passed: twelve runs finished, no
unsampled bins, every overlap above the threshold, no refusal. The nine-window
study before it covered part of the turn only, kept its windows away from the
wrap, and looked entirely reasonable.
"""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

from fastmdxplora.simulation.umbrella import compute_pmf, displacement

KT = 0.008314462618 * 300.0
FORCE = 50.0
BARRIER = 12.0


def _truth(x: np.ndarray) -> np.ndarray:
    """A periodic double well of known height."""
    return BARRIER * (1.0 - np.cos(2.0 * x)) / 2.0


def _study(n_windows: int = 12, seed: int = 0):
    """Windows tiling the whole turn, sampled from the known profile."""
    rng = np.random.default_rng(seed)
    centres = np.linspace(-np.pi, np.pi, n_windows + 1)[:-1]
    grid = np.linspace(-np.pi, np.pi, 4000)
    samples = {}
    for i, centre in enumerate(centres):
        offset = np.remainder(grid - centre + np.pi, 2 * np.pi) - np.pi
        weight = np.exp(-(_truth(grid) + 0.5 * FORCE * offset ** 2) / KT)
        samples[i] = rng.choice(grid, size=3000, p=weight / weight.sum())
    plan = SimpleNamespace(
        windows=[SimpleNamespace(index=i, centre=float(c),
                                 force_constant=FORCE)
                 for i, c in enumerate(centres)],
        collective_variable="torsion", minimum_overlap=0.01,
        minimum_samples=100, equilibration_fraction=0.2,
        equilibration_steps=0)
    return samples, plan


class TestTheDistanceIsMeasuredTheShortWayRound:
    def test_across_the_wrap(self) -> None:
        got = displacement(np.radians([170.0]), np.radians(-180.0), True)
        assert np.degrees(got[0]) == pytest.approx(-10.0, abs=0.1)

    def test_a_straight_coordinate_is_left_alone(self) -> None:
        """A distance has ends. Wrapping one would be nonsense."""
        got = displacement(np.array([3.0]), 0.5, False)
        assert got[0] == pytest.approx(2.5)

    def test_the_energy_that_follows(self) -> None:
        """933 kJ/mol against 0.76 is the whole defect in one number."""
        wrong = 0.5 * FORCE * displacement(
            np.radians([170.0]), np.radians(-180.0), False)[0] ** 2
        right = 0.5 * FORCE * displacement(
            np.radians([170.0]), np.radians(-180.0), True)[0] ** 2
        assert wrong > 900 and right < 1
        assert np.exp((wrong - right) / KT) > 1e100


class TestAKnownPeriodicProfileIsRecovered:
    """The strongest check available: a profile whose height is chosen, and
    windows sampled from it."""

    def test_the_barrier_comes_back(self) -> None:
        out = compute_pmf(*_study(), temperature_K=300.0,
                          bootstrap_resamples=0)
        assert out.get("refused") is None
        energy = np.array([np.nan if v is None else v
                           for v in out["pmf"]["free_energy_kjmol"]])
        finite = energy[np.isfinite(energy)]
        assert (finite.max() - finite.min()) == pytest.approx(BARRIER, abs=1.5)

    def test_the_circle_closes(self) -> None:
        """The two ends are the same point. A chain of window offsets with
        nothing tying them together drifts, and the drift shows up here."""
        out = compute_pmf(*_study(), temperature_K=300.0,
                          bootstrap_resamples=0)
        energy = np.array([np.nan if v is None else v
                           for v in out["pmf"]["free_energy_kjmol"]])
        finite = np.isfinite(energy)
        first, last = energy[finite][0], energy[finite][-1]
        assert abs(first - last) < 2.0, (
            "the profile does not close: the free energy at one end of the "
            "turn differs from the other, which is the same place")

    def test_it_is_not_a_ramp(self) -> None:
        """The failure looked like a plausible file: a monotonic slide with
        no minimum, and every check passing."""
        out = compute_pmf(*_study(), temperature_K=300.0,
                          bootstrap_resamples=0)
        energy = np.array([np.nan if v is None else v
                           for v in out["pmf"]["free_energy_kjmol"]])
        finite = energy[np.isfinite(energy)]
        steps = np.diff(finite)
        assert np.any(steps > 0) and np.any(steps < 0), \
            "a free energy that only ever falls is a drift, not a profile"
