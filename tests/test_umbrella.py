"""A free energy along a coordinate, from equilibrium sampling at each point.

Steered MD drags a system and reports work that depends on how fast it was
dragged. Umbrella sampling holds it at a series of positions, lets each
equilibrate, and recombines the sampling -- so nothing is hurried and nothing
is dissipated.
"""

from __future__ import annotations

import numpy as np
import pytest

from fastmdxplora.simulation.umbrella import KB_KJ

#: A double well with a barrier of exactly 10 kJ/mol, so the recovered value
#: can be checked against a number rather than against plausibility.
def _true_free_energy(x):
    return 10.0 * (np.asarray(x) ** 2 - 1.0) ** 2


def _sample(centre, force, n=6000, seed=0, temperature_K=300.0, width=0.8):
    """Draw from the biased distribution the window would sample."""
    kT = KB_KJ * temperature_K
    rng = np.random.RandomState(seed)
    drawn = []
    while len(drawn) < n:
        x = rng.uniform(centre - width, centre + width, 400)
        energy = _true_free_energy(x) + 0.5 * force * (x - centre) ** 2
        accept = np.exp(-(energy - energy.min()) / kT)
        drawn.extend(x[rng.uniform(size=x.size) < accept])
    return np.array(drawn[:n])


class TestItRecoversAFreeEnergyItWasGiven:
    """The check that matters: sample a landscape whose barrier is known and
    see whether the number comes back."""

    def test_the_barrier_is_recovered(self) -> None:
        from fastmdxplora.simulation.umbrella import compute_pmf, plan_windows

        plan = plan_windows({
            "collective_variable": "distance", "from": -1.6, "to": 1.6,
            "n_windows": 17, "force_constant": 200.0})
        samples = {w.index: _sample(w.centre, w.force_constant, seed=w.index)
                   for w in plan.windows}

        result = compute_pmf(samples, plan, temperature_K=300.0)
        assert result["refused"] is None

        x = np.array(result["pmf"]["coordinate"])
        f = np.array(result["pmf"]["free_energy_kjmol"])
        inside = (x > -1.4) & (x < 1.4) & np.isfinite(f)

        truth = _true_free_energy(x[inside])
        truth = truth - truth.min()
        error = np.sqrt(((f[inside] - truth) ** 2).mean())
        assert error < 1.0, f"recovered PMF is off by {error:.2f} kJ/mol RMS"

    def test_the_barrier_height_is_about_right(self) -> None:
        from fastmdxplora.simulation.umbrella import compute_pmf, plan_windows

        plan = plan_windows({
            "collective_variable": "distance", "from": -1.6, "to": 1.6,
            "n_windows": 17, "force_constant": 200.0})
        result = compute_pmf(
            {w.index: _sample(w.centre, w.force_constant, seed=w.index)
             for w in plan.windows}, plan)

        x = np.array(result["pmf"]["coordinate"])
        f = np.array(result["pmf"]["free_energy_kjmol"])
        at_barrier = f[np.argmin(np.abs(x))]
        assert at_barrier == pytest.approx(10.0, abs=1.5)


class TestWindowsThatDoNotOverlap:
    """Recombination stitches histograms together. Where two neighbours never
    visit the same value there is nothing to stitch, and a curve through the
    gap is interpolation presented as a measurement.
    """

    @staticmethod
    def _separated():
        from fastmdxplora.simulation.umbrella import compute_pmf, plan_windows

        plan = plan_windows({
            "collective_variable": "distance", "from": -1.5, "to": 1.5,
            "n_windows": 4, "force_constant": 3000.0})
        samples = {w.index: _sample(w.centre, w.force_constant, n=4000,
                                    seed=w.index, width=0.5)
                   for w in plan.windows}
        return compute_pmf(samples, plan)

    def test_no_free_energy_is_returned(self) -> None:
        assert self._separated()["pmf"] is None

    def test_the_gap_is_named(self) -> None:
        refusal = self._separated()["refused"]
        assert "do not overlap" in refusal
        assert "windows 0 and 1" in refusal

    def test_and_what_would_close_it(self) -> None:
        """More windows or a softer constant. Not more sampling."""
        refusal = self._separated()["refused"]
        assert "More windows" in refusal
        assert "Sampling for longer will not" in refusal

    def test_the_overlaps_are_reported_either_way(self) -> None:
        """So somebody can see how close the windows came."""
        result = self._separated()
        assert result["overlaps"]
        assert all("overlap" in o for o in result["overlaps"])


class TestPlanningWindows:
    def test_a_range_and_a_count_becomes_windows(self) -> None:
        from fastmdxplora.simulation.umbrella import plan_windows

        plan = plan_windows({"collective_variable": "distance", "from": 0.4,
                             "to": 2.0, "n_windows": 9,
                             "force_constant": 500.0})
        assert len(plan.windows) == 9
        assert plan.windows[0].centre == pytest.approx(0.4)
        assert plan.windows[-1].centre == pytest.approx(2.0)

    def test_explicit_centres_are_taken_as_given(self) -> None:
        from fastmdxplora.simulation.umbrella import plan_windows

        plan = plan_windows({"collective_variable": "distance",
                             "centres": [0.4, 0.5, 0.7, 1.2],
                             "force_constant": 500.0})
        assert [w.centre for w in plan.windows] == [0.4, 0.5, 0.7, 1.2]

    def test_the_force_constant_has_no_default(self) -> None:
        """It sets how far a window wanders and therefore whether neighbours
        overlap. There is no value right for an arbitrary coordinate."""
        from fastmdxplora.simulation.umbrella import plan_windows

        with pytest.raises(ValueError, match="force_constant"):
            plan_windows({"collective_variable": "distance", "from": 0.0,
                          "to": 1.0, "n_windows": 5})

    def test_one_window_is_not_umbrella_sampling(self) -> None:
        from fastmdxplora.simulation.umbrella import plan_windows

        with pytest.raises(ValueError, match="at least two windows"):
            plan_windows({"collective_variable": "distance",
                          "centres": [1.0], "force_constant": 500.0})

    def test_windows_become_a_sweep(self) -> None:
        """One system at many restraint positions is the shape the batch
        machinery already runs, so it schedules them rather than a second
        one."""
        from fastmdxplora.simulation.umbrella import plan_windows, windows_as_sweep

        plan = plan_windows({"collective_variable": "distance", "from": 0.0,
                             "to": 1.0, "n_windows": 5,
                             "force_constant": 500.0})
        runs = windows_as_sweep(plan)
        assert len(runs) == 5
        assert runs[0]["id"] == "window_00"
        assert runs[-1]["umbrella_centre"] == pytest.approx(1.0)


class TestOverlapItself:
    def test_identical_sampling_overlaps_completely(self) -> None:
        from fastmdxplora.simulation.umbrella import overlap_between

        x = np.random.RandomState(0).normal(size=2000)
        assert overlap_between(x, x) == pytest.approx(1.0, abs=0.05)

    def test_separated_sampling_does_not_overlap(self) -> None:
        from fastmdxplora.simulation.umbrella import overlap_between

        rng = np.random.RandomState(0)
        assert overlap_between(rng.normal(0, 0.1, 2000),
                               rng.normal(5, 0.1, 2000)) < 0.01

    def test_empty_sampling_overlaps_nothing(self) -> None:
        from fastmdxplora.simulation.umbrella import overlap_between

        assert overlap_between(np.array([]), np.array([1.0, 2.0])) == 0.0
