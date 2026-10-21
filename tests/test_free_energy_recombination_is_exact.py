"""A free energy recombined from windows whose answer is known exactly.

Windows drawn from the exact biased distributions of a known potential of
mean force, so every number can be checked against a closed form rather than
against itself. On that test a binding free energy was refused two times in
three: the curve ran past the outermost windows into bins one window's tail
reached, solved as though the density were flat across a bin, and the bulk
check read exactly those bins. Where it was not refused it was 0.2 kJ/mol off.
"""

from __future__ import annotations

import json
import logging

import numpy as np
import pytest

from fastmdxplora.simulation.binding import (
    KB_KJMOL, STANDARD_VOLUME_NM3, binding_free_energy)
from fastmdxplora.simulation.umbrella import (
    UmbrellaPlan, Window, compute_pmf, describe_pmf)

TEMPERATURE = 300.0
KT = KB_KJMOL * TEMPERATURE
GRID = np.linspace(0.2, 2.2, 200001)
#: The true potential of mean force: a 20 kJ/mol well at 0.5 nm, flat beyond.
W_TRUE = -20.0 * np.exp(-((GRID - 0.5) / 0.1) ** 2)
#: The free energy of the distance's distribution, the shell's room in it.
A_TRUE = W_TRUE - 2.0 * KT * np.log(GRID)
CUTOFF = 0.8


def _exact_binding() -> float:
    inside = GRID <= CUTOFF
    constant = 4.0 * np.pi * np.trapezoid(
        GRID[inside] ** 2 * np.exp(-W_TRUE[inside] / KT), GRID[inside])
    return float(-KT * np.log(constant / STANDARD_VOLUME_NM3))


def _plan(variable="ligand_distance", centres=None, force=1000.0):
    centres = np.linspace(0.3, 2.0, 36) if centres is None else centres
    return UmbrellaPlan(
        windows=tuple(Window(index=i, centre=float(c), force_constant=force)
                      for i, c in enumerate(centres)),
        collective_variable=variable, equilibration_fraction=0.2,
        minimum_overlap=0.0, minimum_samples=10, cone=None,
        bootstrap_resamples=0)


def _samples(plan, count, seed):
    rng = np.random.default_rng(seed)
    drawn = {}
    for window in plan.windows:
        weight = np.exp(-(A_TRUE + 0.5 * window.force_constant
                          * (GRID - window.centre) ** 2) / KT)
        drawn[window.index] = rng.choice(GRID, size=count,
                                         p=weight / weight.sum())
    return drawn


@pytest.fixture(scope="module")
def recombined():
    logging.disable(logging.WARNING)
    try:
        plan = _plan()
        return [compute_pmf(_samples(plan, 5000, 100 + seed), plan,
                            temperature_K=TEMPERATURE, bootstrap_resamples=0)
                for seed in range(6)]
    finally:
        logging.disable(logging.NOTSET)


class TestTheBindingFreeEnergyIsTheClosedForm:
    def test_it_is_not_refused_and_lands_on_it(self, recombined) -> None:
        exact = _exact_binding()
        found = []
        for payload in recombined:
            got = binding_free_energy(
                payload["pmf"]["coordinate"],
                payload["pmf"]["free_energy_kjmol"],
                temperature_K=TEMPERATURE, bound_cutoff_nm=CUTOFF)
            assert got["refused"] is None, got["refused"]
            found.append(got["delta_g_kjmol"])
        # Six draws of 5000 samples a window: the mean within a tenth of a
        # kJ/mol of the closed form (it was 0.24 off before, when not refused).
        assert np.mean(found) == pytest.approx(exact, abs=0.1)

    def test_the_curve_stops_at_the_last_window(self, recombined) -> None:
        # Outward, at the last window's centre, where the bulk is read; inward,
        # where the data begin, so the well's inner side is kept.
        x = np.asarray(recombined[0]["pmf"]["coordinate"])
        assert 1.98 < x.max() <= 2.0
        assert x.min() < 0.3
        assert recombined[0]["unsampled_bins"] == 0


def _walled(r):
    """The well with a wall the ligand cannot enter below about 0.36 nm."""
    return (-20.0 * np.exp(-((r - 0.5) / 0.1) ** 2)
            + 10.0 * (0.36 / r) ** 12)


class TestWindowsPlacedWhereTheLigandCannotGo:
    """Windows begun inside the wall, or at the well itself: the curve once
    started at the first window's centre, opened with bins no sample could
    reach, and refused the binding free energy as a gap more sampling would
    close."""

    @pytest.mark.parametrize("first", [0.25, 0.5])
    def test_it_is_not_refused_and_lands_on_it(self, first) -> None:
        w_true = _walled(GRID)
        a_true = w_true - 2.0 * KT * np.log(GRID)
        inside = GRID <= CUTOFF
        exact = float(-KT * np.log(4.0 * np.pi * np.trapezoid(
            GRID[inside] ** 2 * np.exp(-w_true[inside] / KT), GRID[inside])
            / STANDARD_VOLUME_NM3))
        plan = _plan(centres=np.linspace(first, 2.0, 36))
        rng = np.random.default_rng(7)
        drawn = {}
        for window in plan.windows:
            weight = np.exp(-(a_true + 0.5 * window.force_constant
                              * (GRID - window.centre) ** 2) / KT)
            drawn[window.index] = rng.choice(GRID, size=5000,
                                             p=weight / weight.sum())
        logging.disable(logging.WARNING)
        try:
            payload = compute_pmf(drawn, plan, temperature_K=TEMPERATURE,
                                  bootstrap_resamples=0)
        finally:
            logging.disable(logging.NOTSET)
        got = binding_free_energy(payload["pmf"]["coordinate"],
                                  payload["pmf"]["free_energy_kjmol"],
                                  temperature_K=TEMPERATURE,
                                  bound_cutoff_nm=CUTOFF)
        assert got["refused"] is None, got["refused"]
        assert got["delta_g_kjmol"] == pytest.approx(exact, abs=0.6)


class TestTheCurveAlongADistanceIsSaidAsWhatItIs:
    def test_w_of_r_is_flat_in_bulk_and_a_is_not(self, recombined) -> None:
        pmf = recombined[0]["pmf"]
        x = np.asarray(pmf["coordinate"])
        w = np.asarray(pmf["potential_of_mean_force_kjmol"], dtype=float)
        a = np.asarray(pmf["free_energy_kjmol"], dtype=float)
        bulk = x > 1.2
        assert np.ptp(w[bulk]) < 1.0
        assert np.ptp(a[bulk]) > 2.0  # 2kT ln(2.0 / 1.2) = 2.5
        assert pmf["form"] == "distribution"

    def test_a_torsion_has_no_shell(self) -> None:
        centres = np.linspace(-2.5, 2.5, 12)
        plan = _plan("torsion", centres=centres, force=50.0)
        rng = np.random.default_rng(0)
        drawn = {w.index: w.centre + rng.normal(scale=0.2, size=400)
                 for w in plan.windows}
        payload = compute_pmf(drawn, plan, temperature_K=TEMPERATURE,
                              bootstrap_resamples=0)
        assert payload["pmf"]["potential_of_mean_force_kjmol"] is None
        assert payload["pmf"]["form"] == "coordinate"


def _a_curve_with_a_hole(at_nm):
    x = np.linspace(0.3, 2.0, 60)
    a = (-20.0 * np.exp(-((x - 0.5) / 0.1) ** 2) - 2.0 * KT * np.log(x))
    energy = [float(v) for v in a - a.min()]
    energy[int(np.argmin(np.abs(x - at_nm)))] = None
    return x, energy


class TestAnEmptyBinIsSaidWhereItIs:
    def test_in_the_bound_range(self) -> None:
        x, energy = _a_curve_with_a_hole(0.45)
        got = binding_free_energy(x, energy, temperature_K=TEMPERATURE)
        assert got["delta_g_kjmol"] is None
        assert "inside the bound range" in got["refused"]
        assert "no population above the bulk reference" not in got["refused"]

    def test_in_the_bulk_range(self) -> None:
        x, energy = _a_curve_with_a_hole(1.9)
        got = binding_free_energy(x, energy, temperature_K=TEMPERATURE)
        assert got["delta_g_kjmol"] is None
        assert "bulk reference is read from" in got["refused"]
        # Refused, so nothing goes on to check the reference.
        assert "bulk_from_nm" not in got

    def test_at_the_wall_it_holds_nobody(self) -> None:
        x = np.linspace(0.3, 2.0, 60)
        a = _walled(x) - 2.0 * KT * np.log(x)
        energy = [float(v) for v in a - a.min()]
        energy[1] = None
        got = binding_free_energy(x, energy, temperature_K=TEMPERATURE,
                                  bound_cutoff_nm=CUTOFF)
        assert got["refused"] is None

    def test_between_them_it_costs_nothing(self) -> None:
        x, energy = _a_curve_with_a_hole(1.1)
        got = binding_free_energy(x, energy, temperature_K=TEMPERATURE,
                                  bound_cutoff_nm=CUTOFF)
        assert got["refused"] is None
        assert np.isfinite(got["delta_g_kjmol"])


class TestTheAnalysisDrawsW:
    def test_an_older_record_gets_its_w_worked_out(self, tmp_path) -> None:
        from fastmdxplora.analysis.pmf import PMF

        x = np.linspace(0.3, 2.0, 60)
        a = np.interp(x, GRID, A_TRUE)
        (tmp_path / "pmf.json").write_text(json.dumps({
            "pmf": {"coordinate": list(x), "free_energy_kjmol": list(a - a.min())},
            "plan": {"collective_variable": "ligand_distance"},
            "temperature_K": TEMPERATURE, "covered": [0.3, 2.0],
            "summary": describe_pmf(x, list(a - a.min()), (0.3, 2.0)),
            "refused": None, "unsampled_bins": 0, "n_windows": 36,
        }), encoding="utf-8")
        analysis = PMF(output_dir=tmp_path / "analysis")
        result = analysis.compute(None)
        w = result["potential_of_mean_force_kjmol"]
        assert w is not None and np.ptp(w[x > 1.2]) < 0.2
        # Its summary read again off W: a plain well has no barrier, and its
        # highest point, the plateau, is the well's depth, not the 22.2
        # kJ/mol of A's wall at the first window.
        assert result["summary"]["barrier_kjmol"] is None
        assert result["summary"]["highest_kjmol"] == pytest.approx(20.0, abs=0.3)
        assert analysis.default_ylabel().startswith("W(r)")
        written = analysis.save_data(result, tmp_path / "pmf.dat").read_text()
        assert "potential_of_mean_force_kjmol" in written.splitlines()[0]


class TestABarrierIsAPeak:
    """The barrier was the curve's highest point: along a distance the wall
    at the first window, on a full turn tiled the usual way nothing at
    +-pi, and in a plain well under soft windows a bin few samples reached."""

    def test_a_plain_well_recombined_has_none(self, recombined) -> None:
        for payload in recombined:
            assert payload["summary"]["barrier_kjmol"] is None
            assert 19.0 < payload["summary"]["highest_kjmol"] < 20.6

    def test_a_hump_inside_is_one(self) -> None:
        x = np.linspace(0.0, 1.0, 41)
        e = 30.0 * (1.0 - x) + 8.0 * np.exp(-((x - 0.6) / 0.05) ** 2)
        got = describe_pmf(x, list(e))
        assert got["barrier_at"] == pytest.approx(0.6, abs=0.03)
        assert got["highest_is_an_end"] is True

    def test_a_spike_within_its_error_is_not(self) -> None:
        x = np.linspace(0.3, 2.0, 60)
        e = np.where(x > 0.7, 20.0, 0.0) + 0.0 * x
        e[20] += 4.0
        quiet = describe_pmf(x, list(e))
        assert quiet["barrier_kjmol"] is not None
        errors = np.full(x.size, 0.2)
        errors[20] = 3.0
        assert describe_pmf(x, list(e), errors=list(errors))["barrier_kjmol"] is None

    @pytest.mark.parametrize("first", [-np.pi, -np.pi + np.pi / 24])
    def test_one_at_the_turn_s_ends_is_found(self, first) -> None:
        # 24 windows round the turn, none of them on +pi: the slice they
        # cover falls short of a full turn, and a peak rule judging closure
        # on that slice lost the barrier at +-pi.
        centres = first + np.arange(24) * (2.0 * np.pi / 24)
        plan = _plan("torsion", centres=centres, force=200.0)
        grid = np.linspace(-np.pi, np.pi, 20001)
        truth = 7.5 * (1.0 - np.cos(grid))  # the barrier at +-pi
        rng = np.random.default_rng(11)
        drawn = {}
        for window in plan.windows:
            d = np.remainder(grid - window.centre + np.pi, 2.0 * np.pi) - np.pi
            weight = np.exp(-(truth + 0.5 * window.force_constant * d ** 2) / KT)
            drawn[window.index] = rng.choice(grid, size=3000, p=weight / weight.sum())
        logging.disable(logging.WARNING)
        try:
            payload = compute_pmf(drawn, plan, temperature_K=TEMPERATURE,
                                  bootstrap_resamples=0)
        finally:
            logging.disable(logging.NOTSET)
        summary = payload["summary"]
        assert summary["barrier_kjmol"] == pytest.approx(15.0, abs=1.0)
        assert abs(abs(summary["barrier_at"]) - np.pi) < 0.3


def _torsion_windows(centres, force, samples, truth_of, seed):
    plan = _plan("torsion", centres=np.asarray(centres), force=force)
    grid = np.linspace(-np.pi, np.pi, 20001)
    truth = truth_of(grid)
    rng = np.random.default_rng(seed)
    drawn = {}
    for window in plan.windows:
        d = np.remainder(grid - window.centre + np.pi, 2.0 * np.pi) - np.pi
        weight = np.exp(-(truth + 0.5 * window.force_constant * d ** 2) / KT)
        drawn[window.index] = rng.choice(grid, size=samples, p=weight / weight.sum())
    logging.disable(logging.WARNING)
    try:
        return compute_pmf(drawn, plan, temperature_K=TEMPERATURE,
                           bootstrap_resamples=0)
    finally:
        logging.disable(logging.NOTSET)


class TestAPartialArcIsNotATurn:
    def test_tails_that_wrap_do_not_close_it(self) -> None:
        # Nine windows from -30 to 165 degrees: their tails reached round to
        # +-pi, and judged from the samples the arc read as a full turn, its
        # minimum moved into the arc no window held.
        payload = _torsion_windows(np.radians(np.linspace(-30, 165, 9)), 200.0,
                                   3000, lambda x: 9.0 * (1.0 - np.cos(3 * x)) / 2.0,
                                   seed=5)
        low, high = np.degrees(payload["summary"]["covered"])
        assert -31 < low and high < 166


class TestAnOlderRecordIsReadAgain:
    def test_its_barrier_is_a_peak(self, tmp_path) -> None:
        from fastmdxplora.analysis.pmf import PMF

        x = np.linspace(-np.pi, np.pi, 60)
        a = 7.5 * (1.0 - np.cos(x))
        (tmp_path / "pmf.json").write_text(json.dumps({
            "pmf": {"coordinate": list(x), "free_energy_kjmol": list(a),
                    "uncertainty": {"standard_error": [0.2] * 60}},
            "plan": {"collective_variable": "torsion", "n_windows": 24},
            "temperature_K": TEMPERATURE,
            "covered": [-np.pi, np.pi - 2.0 * np.pi / 24],
            # Written by the earlier rule: the highest point within covered.
            "summary": {"barrier_kjmol": 14.8, "covered": [-3.1, 2.9],
                        "minima": []},
            "refused": None, "unsampled_bins": 0, "n_windows": 24,
        }), encoding="utf-8")
        result = PMF(output_dir=tmp_path / "analysis").compute(None)
        assert result["summary"]["barrier_kjmol"] == pytest.approx(15.0, abs=0.1)
        assert abs(abs(result["summary"]["barrier_at"]) - np.pi) < 0.1


class TestATurnKeepsItsEndsAsWells:
    def test_the_trans_well_at_the_ends_is_a_minimum(self) -> None:
        x = np.linspace(-np.pi, np.pi, 72, endpoint=False) + np.pi / 72
        e = 6.0 * (1.0 + np.cos(x)) + 3.0 * (1.0 + np.cos(3.0 * x))
        got = describe_pmf(x, list(e), periodic=True, closed=True)
        deepest = got["minima"][0]
        assert deepest["free_energy_kjmol"] == pytest.approx(0.0, abs=0.1)
        assert abs(abs(deepest["coordinate"]) - np.pi) < 0.1
        assert got["barrier_kjmol"] == pytest.approx(18.0, abs=0.5)


class TestAnArcAcrossTheSeamIsAnArc:
    def test_the_widest_gap_decides(self) -> None:
        from fastmdxplora.simulation.umbrella import turn_is_tiled

        # 120 to 240 degrees, listed across +-180: not a turn.
        across = list(np.radians(np.r_[np.arange(120, 181, 15), np.arange(-165, -119, 15)]))
        assert not turn_is_tiled(min(across), max(across), len(across), across)
        # 24 windows round the turn: a turn.
        whole = list(np.radians(np.arange(-180, 180, 15)))
        assert turn_is_tiled(min(whole), max(whole), len(whole), whole)
        # A ring with a 90 degree hole about 0: not a turn.
        holed = [c for c in whole if abs(np.degrees(c)) >= 45]
        assert not turn_is_tiled(min(holed), max(holed), len(holed), holed)
