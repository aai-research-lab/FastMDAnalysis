"""A metadynamics run produced no result at all.

It wrote HILLS and COLVAR and stopped: the module said "a run that has not
converged has no free energy" while offering no free energy either way, so
somebody wanting a surface ran ``plumed sum_hills`` themselves and got one
with nothing attached to it. Umbrella sampling went from a config block to a
curve or a refusal; metadynamics went to a pair of PLUMED files.
"""

from __future__ import annotations

import numpy as np
import pytest

from fastmdxplora.simulation.metad_surface import (
    SETTLED_DRIFT_KJMOL,
    Hills,
    compute_surface,
    read_hills,
    recrossings,
    surface_from_hills,
)

GAMMA = 10.0
SIGMA = 0.05


def _grid():
    return np.linspace(0.0, 1.0, 400)


def _double_well(grid):
    """Two basins with a 20 kJ/mol barrier between them."""
    surface = 20.0 * (1.0 - np.cos(2.0 * np.pi * grid)) / 2.0
    return surface - surface.min()


def _hills_reaching(free_energy, grid, *, bias_factor=GAMMA, centres=400):
    """Hills whose sum is the well-tempered bias for the given free energy.

    Heights by quadrature -- each hill carries the bias at its centre, times
    the spacing over the Gaussian's normalisation -- so they are positive and
    track the bias, which is what deposition actually produces. A first
    version solved for the heights by least squares: the sum was exact, and
    the individual heights oscillated to plus and minus a hundred thousand
    and cancelled, so any partial sum was nonsense.

    Constructed rather than simulated, because what is checked here is the
    arithmetic taking a bias to a free energy. A homemade sampler would be
    checking the sampler: the first attempt read 27.7 kJ/mol for a 20 kJ/mol
    barrier because it wrapped the coordinate and not the bias grid. The
    physics is checked against PLUMED on alanine dipeptide, not here.
    """
    scale = (1.0 - 1.0 / bias_factor) if bias_factor > 1.0 else 1.0
    # Deposition fills wells, so the bias is high where the free energy is
    # low. The constant drops out when the surface is shifted.
    bias = scale * (free_energy.max() - free_energy)

    # Placed past both ends, because a hill at the boundary has no neighbours
    # outside it and the sum falls away there: confined to the grid, a 20
    # kJ/mol barrier came back as 16.9.
    margin = 4.0 * SIGMA
    places = np.linspace(grid[0] - margin, grid[-1] + margin, centres)
    spacing = places[1] - places[0]
    heights = (np.interp(places, grid, bias) * spacing
               / (SIGMA * np.sqrt(2.0 * np.pi)))
    return Hills(
        time_ps=np.arange(centres, dtype=float),
        centre=places,
        sigma=np.full(centres, SIGMA),
        height=heights,
        bias_factor=bias_factor,
    )


def _written(tmp_path, hills, name="HILLS"):
    path = tmp_path / name
    path.write_text(
        "#! FIELDS time cv sigma_cv height biasf\n"
        + "\n".join(
            f"{t:.3f} {c:.6f} {s:.6f} {h:.6f} {hills.bias_factor:.1f}"
            for t, c, s, h in zip(hills.time_ps, hills.centre,
                                  hills.sigma, hills.height))
        + "\n",
        encoding="utf-8")
    return path


class TestTheBiasBecomesAFreeEnergy:
    """Well-tempered deposition converges on -(1 - 1/γ)F, so the free energy
    is the bias scaled by γ/(γ - 1) and negated."""

    def test_a_known_barrier_comes_back(self) -> None:
        """To within the smoothing the hills impose -- see below."""
        grid = _grid()
        expected = _double_well(grid)
        hills = _hills_reaching(expected, grid)

        recovered = surface_from_hills(hills, grid)

        assert recovered.max() == pytest.approx(expected.max(), rel=0.05)

    def test_the_hills_smooth_the_surface_they_build(self) -> None:
        """A sum of Gaussians of width sigma is the surface convolved with
        one, so a feature narrower than sigma comes back shallower than it is.
        A period-1 cosine is damped by exp(-2 pi^2 sigma^2), which for sigma
        0.05 is 0.95 -- and that is the whole reason sigma has to be smaller
        than the features being resolved, rather than a number to leave at a
        default.

        Asserted rather than tolerated, because it looked like an error of
        4% until the arithmetic said otherwise.
        """
        grid = _grid()
        expected = _double_well(grid)
        damping = float(np.exp(-2.0 * np.pi ** 2 * SIGMA ** 2))

        recovered = surface_from_hills(_hills_reaching(expected, grid), grid)

        assert recovered.max() == pytest.approx(
            expected.max() * damping, rel=0.02)

    def test_a_wider_hill_smooths_more(self) -> None:
        """Which is the same statement, made where it can be seen."""
        grid = _grid()
        expected = _double_well(grid)

        global SIGMA
        narrow = surface_from_hills(_hills_reaching(expected, grid), grid).max()
        was, SIGMA = SIGMA, 0.12
        try:
            wide = surface_from_hills(_hills_reaching(expected, grid), grid).max()
        finally:
            SIGMA = was

        assert wide < narrow < expected.max()

    def test_without_tempering_the_bias_is_the_free_energy(self) -> None:
        """A bias factor of one means no tempering, and the scaling would
        divide by zero."""
        grid = _grid()
        expected = _double_well(grid)
        hills = _hills_reaching(expected, grid, bias_factor=1.0)

        recovered = surface_from_hills(hills, grid)
        assert recovered.max() == pytest.approx(expected.max(), rel=0.05)

    def test_the_lowest_point_is_zero(self) -> None:
        """Only differences mean anything: where the sum started is an
        accident of the run."""
        grid = _grid()
        hills = _hills_reaching(_double_well(grid), grid)
        assert surface_from_hills(hills, grid).min() == pytest.approx(0.0)

    def test_it_can_be_built_from_part_of_the_run(self) -> None:
        """Which is what makes the drift check possible."""
        grid = _grid()
        hills = _hills_reaching(_double_well(grid), grid)

        whole = surface_from_hills(hills, grid)
        partial = surface_from_hills(hills, grid, upto=len(hills) // 2)

        # Half the hills describe half the coordinate, so the two disagree --
        # which is the quantity the drift check measures.
        assert np.max(np.abs(whole - partial)) > 1.0


class TestCountingWhenTheSystemCameBack:
    """A barrier crossed once has been observed once, and its height rests on
    that single event."""

    def test_a_round_trip_counts_twice(self) -> None:
        values = np.array([0.0, 0.5, 1.0, 0.5, 0.0])
        assert recrossings(values, low=0.25, high=0.75) == 2

    def test_jitter_at_a_threshold_counts_for_nothing(self) -> None:
        """Without hysteresis a coordinate rattling at the top of a barrier
        accumulates crossings it never made."""
        values = np.array([0.74, 0.76, 0.74, 0.76, 0.74, 0.76] * 20)
        assert recrossings(values, low=0.25, high=0.75) == 0

    def test_a_coordinate_that_never_returns_counts_once(self) -> None:
        values = np.concatenate([np.zeros(100), np.ones(100)])
        assert recrossings(values, low=0.25, high=0.75) == 1

    def test_nothing_sampled_is_no_crossings(self) -> None:
        assert recrossings(np.array([]), low=0.25, high=0.75) == 0


class TestARunThatDoesNotSupportASurface:
    """Three conditions, each answerable from the files the run already
    writes, and each refused by name."""

    def _converged(self, tmp_path):
        grid = _grid()
        hills = _hills_reaching(_double_well(grid), grid)
        # Heights that have decayed deep into the tail, as well-tempered
        # deposition does once a basin is full.
        hills = Hills(hills.time_ps, hills.centre, hills.sigma,
                      hills.height, hills.bias_factor)
        return _written(tmp_path, hills), grid

    def test_hills_still_arriving_at_full_height_are_refused(
        self, tmp_path
    ) -> None:
        grid = _grid()
        flat = _hills_reaching(_double_well(grid), grid)
        # Every hill the same height: nothing has decayed, so the bias is
        # still filling the landscape.
        same = Hills(flat.time_ps, flat.centre, flat.sigma,
                     np.full(flat.centre.size, 1.2), flat.bias_factor)
        crossing = np.tile([0.0, 1.0], 20)

        result = compute_surface(_written(tmp_path, same), crossing)

        assert result["surface"] is None
        assert "still arriving" in result["refused"]

    def test_a_surface_still_moving_is_refused(self, tmp_path) -> None:
        """Built from three quarters of the hills and from all of them, a
        converged surface gives the same answer."""
        grid = _grid()
        hills = _hills_reaching(_double_well(grid), grid)
        result = compute_surface(_written(tmp_path, hills), np.tile([0., 1.], 20))

        # These hills are deposited left to right rather than by sampling, so
        # the last quarter carries a whole basin: the surface moves a great
        # deal and the check says so.
        assert result["surface"] is None
        assert "has not stopped changing" in result["refused"]
        assert result["evidence"]["drift_kjmol"] > SETTLED_DRIFT_KJMOL

    def test_too_few_crossings_are_refused_by_number(self, tmp_path) -> None:
        grid = _grid()
        hills = _hills_reaching(_double_well(grid), grid)
        one_way = np.concatenate([np.zeros(50), np.ones(50)])

        refused = compute_surface(_written(tmp_path, hills), one_way)["refused"]

        assert "crossed between the ends of the range 1 time(s)" in refused
        assert "rests on that single event" in refused

    def test_no_trajectory_is_said_rather_than_skipped(self, tmp_path) -> None:
        """A surface reported without knowing whether the system ever came
        back is the claim this exists to avoid."""
        grid = _grid()
        hills = _hills_reaching(_double_well(grid), grid)

        refused = compute_surface(_written(tmp_path, hills), None)["refused"]
        assert "was not available" in refused

    def test_a_coordinate_that_never_moved_is_refused(self, tmp_path) -> None:
        still = Hills(np.arange(50.0), np.full(50, 0.3), np.full(50, SIGMA),
                      np.full(50, 1.0), GAMMA)
        result = compute_surface(_written(tmp_path, still), np.full(50, 0.3))

        assert result["surface"] is None
        assert "did not move" in result["refused"]

    def test_the_evidence_is_returned_with_the_refusal(self, tmp_path) -> None:
        """So a study can be read from its files rather than from a sentence,
        and a borderline run can be judged rather than only rejected."""
        grid = _grid()
        hills = _hills_reaching(_double_well(grid), grid)
        result = compute_surface(_written(tmp_path, hills), np.tile([0., 1.], 20))

        assert set(result["evidence"]) == {
            "hills", "bias_factor", "first_hill_height_kjmol",
            "last_hill_height_kjmol", "drift_kjmol", "recrossings",
            "barrier_kjmol"}

    def test_every_refusal_says_that_longer_is_the_remedy(self, tmp_path) -> None:
        """Because it is, for all three -- unlike an umbrella gap, where
        sampling for longer will not close it."""
        grid = _grid()
        hills = _hills_reaching(_double_well(grid), grid)
        refused = compute_surface(_written(tmp_path, hills), None)["refused"]
        assert "Sampling for longer is the remedy" in refused


class TestReadingWhatPlumedWrote:
    def test_comments_and_blank_lines_are_skipped(self, tmp_path) -> None:
        path = tmp_path / "HILLS"
        path.write_text(
            "#! FIELDS time cv sigma_cv height biasf\n"
            "#! SET multivariate false\n"
            "\n"
            "0.000 0.100000 0.050000 1.200000 10.0\n"
            "1.000 0.200000 0.050000 1.100000 10.0\n", encoding="utf-8")

        hills = read_hills(path)
        assert len(hills) == 2
        assert hills.bias_factor == 10.0
        assert hills.centre[1] == pytest.approx(0.2)

    def test_a_file_with_no_hills_is_refused_not_summed(self, tmp_path) -> None:
        path = tmp_path / "HILLS"
        path.write_text("#! FIELDS time cv sigma_cv height biasf\n",
                        encoding="utf-8")

        with pytest.raises(ValueError, match="deposited nothing"):
            read_hills(path)


class TestTheRunActuallyProducesIt:
    """The umbrella pieces were once committed unreachable -- the expansion
    called by nothing, then the recombination called by nothing -- so a config
    would have been accepted, ignored and run once. A surface nothing builds
    is the same defect.
    """

    def test_the_simulation_phase_builds_one(self) -> None:
        import inspect

        from fastmdxplora.simulation import pipeline

        source = inspect.getsource(pipeline.run)
        assert '_write_metadynamics_surface(output_dir' in source
        assert 'params.get("metadynamics")' in source

    def test_and_only_for_a_metadynamics_run(self) -> None:
        """A steered or umbrella run writes a COLVAR too, and summing hills
        that are not there is not a question worth asking."""
        import inspect

        from fastmdxplora.simulation import pipeline

        source = inspect.getsource(pipeline.run)
        guarded = source[source.index('if params.get("metadynamics")'):]
        assert "_write_metadynamics_surface" in guarded.split("\n\n")[0]

    def test_a_refusal_is_written_down_too(self, tmp_path) -> None:
        """A refusal that leaves no trace is a run that looks as though the
        question was never asked."""
        import json

        from fastmdxplora.simulation.pipeline import _write_metadynamics_surface

        grid = _grid()
        hills = _hills_reaching(_double_well(grid), grid)
        _written(tmp_path, hills)
        (tmp_path / "COLVAR").write_text(
            "#! FIELDS time cv\n" + "\n".join(f"{i} 0.3" for i in range(50)),
            encoding="utf-8")

        name = _write_metadynamics_surface(tmp_path, None)

        assert name == "metadynamics_surface.json"
        record = json.loads((tmp_path / name).read_text(encoding="utf-8"))
        assert record["refused"]
        assert record["free_energy_kjmol"] is None
        assert record["evidence"]["recrossings"] == 0

    def test_nothing_is_written_where_there_were_no_hills(self, tmp_path) -> None:
        from fastmdxplora.simulation.pipeline import _write_metadynamics_surface

        assert _write_metadynamics_surface(tmp_path, None) is None


class TestARefusalIsPrintedWhole:
    """The message explaining why there is no surface was cut at 160
    characters, mid-sentence: "...rather than having flattened it -- the
    surface". What went missing was the clause saying the output is still a
    usable snapshot of the filling."""

    def test_nothing_clips_it(self) -> None:
        import inspect

        from fastmdxplora.simulation import pipeline

        source = inspect.getsource(pipeline)
        refusal = source.split('"No free energy surface: "', 1)[1][:300]
        assert "[:160]" not in refusal
        assert "[:" not in refusal.split("presenter.step")[0]

    def test_the_reason_says_the_snapshot_is_still_there(self) -> None:
        """The clause the cut removed."""
        import inspect

        from fastmdxplora.simulation import metad_surface

        source = inspect.getsource(metad_surface)
        assert "is a snapshot of that filling" in source
