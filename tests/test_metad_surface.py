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

        # The snapshot is returned, marked provisional. It used to be None,
        # which contradicted the refusal's own words -- "the surface below is
        # a snapshot of that filling" -- and left a run whose wells had
        # converged with nothing to plot. "Not converged" and "not available"
        # are different claims.
        assert result["surface"] is not None
        assert result["provisional"] is True
        assert "still arriving" in result["refused"]

    def test_a_surface_still_moving_is_refused(self, tmp_path) -> None:
        """Built from three quarters of the hills and from all of them, a
        converged surface gives the same answer."""
        grid = _grid()
        hills = _hills_reaching(_double_well(grid), grid)
        result = compute_surface(_written(tmp_path, hills), np.tile([0., 1.], 20))

        # These hills are deposited left to right rather than by sampling, so
        # the last quarter carries a whole basin: the surface moves a great
        # deal and the check says so -- while still handing back what it has.
        assert result["refused"]
        assert result["surface"] is not None
        assert result["provisional"] is True
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
            # Where the two states are, because the count is of travel
            # between them and a reader cannot check it otherwise.
            "basins",
            # What "recrossings" counts, because the word is ambiguous and
            # the two readings differed by 60% on a real run.
            "recrossings_definition",
            # Where drift was judged, and what it would have read over the
            # whole grid: 2.3 against 5.4 on a real torsion, the difference
            # being one point on top of a 65 kJ/mol barrier.
            "drift_ceiling_kjmol", "drift_over_the_whole_grid_kjmol",
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
        # Written, not withheld: the refusal says the surface below is a
        # snapshot of the filling, and a record that carries the sentence
        # without the numbers leaves a reader with nothing to look at.
        assert record["free_energy_kjmol"] is not None
        assert record["provisional"] is True
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


class TestTheRecrossingCountSaysWhatItCounts:
    """"Recrossings" is a word whose definition changes the number. A run
    recorded as 61 gave 97 by counting sign changes of the raw coordinate --
    60% apart -- and anyone comparing their own count against the record, not
    knowing which was meant, would draw the wrong conclusion about their
    sampling."""

    def test_jitter_is_not_a_traversal(self) -> None:
        """The distinction the hysteresis band exists for: a coordinate
        rattling at a threshold scores hundreds by sign changes and none
        here."""
        import numpy as np

        from fastmdxplora.simulation.metad_surface import recrossings

        jitter = np.random.default_rng(0).normal(0.0, 0.05, 400)
        assert recrossings(jitter, low=-1.5, high=1.5) == 0
        # What a naive count would have said.
        assert int(np.sum(np.diff(np.sign(jitter)) != 0)) > 100

    def test_a_genuine_traversal_counts(self) -> None:
        import numpy as np

        from fastmdxplora.simulation.metad_surface import recrossings

        back_and_forth = np.tile(
            np.concatenate([np.full(20, -3.0), np.full(20, 3.0)]), 10)
        assert recrossings(back_and_forth, low=-1.5, high=1.5) == 19

    def test_the_record_carries_the_definition(self) -> None:
        import inspect

        from fastmdxplora.simulation import metad_surface

        source = inspect.getsource(metad_surface)
        assert '"recrossings_definition"' in source
        # It quotes the band, so the number can be reproduced.
        assert "low_edge" in source and "high_edge" in source


class TestDriftIsJudgedWhereTheSurfaceMeansSomething:
    """Comparing every point on the grid makes the number the worst point
    rather than the typical one, and the worst point is always the top of the
    highest barrier -- estimated from a handful of visits out of thousands of
    hills.

    A tripeptide's psi torsion settled to 1.2 kJ/mol within 10 of its
    minimum and 2.3 within 20, while the whole-grid figure read 5.4 from a
    single point atop a 65 kJ/mol barrier. On that measure no steep
    coordinate can pass however well its wells are resolved: the test could
    not say yes to a correct answer."""

    def test_the_ceiling_is_recorded_with_the_number(self) -> None:
        """A drift measured over part of the range means nothing without
        saying which part."""
        grid = _grid()
        hills = _hills_reaching(_double_well(grid), grid)
        import tempfile
        from pathlib import Path

        result = compute_surface(
            _written(Path(tempfile.mkdtemp()), hills), np.tile([0., 1.], 20))
        evidence = result["evidence"]

        assert "drift_ceiling_kjmol" in evidence
        assert "drift_over_the_whole_grid_kjmol" in evidence
        assert evidence["drift_kjmol"] <= evidence["drift_over_the_whole_grid_kjmol"]

    def test_the_ceiling_is_a_few_multiples_of_rt(self) -> None:
        """A state rarer than this is not one a simulation is measuring."""
        from fastmdxplora.simulation.metad_surface import DRIFT_CEILING_KJMOL

        rt_300k = 2.4789
        assert 4 * rt_300k < DRIFT_CEILING_KJMOL < 15 * rt_300k

    def test_a_barrier_top_does_not_decide_a_converged_run(self) -> None:
        """Wells stable, one high point moving: the verdict should follow the
        wells."""
        import tempfile
        from pathlib import Path

        from fastmdxplora.simulation.metad_surface import (
            DRIFT_CEILING_KJMOL,
            SETTLED_DRIFT_KJMOL,
        )

        # Built directly: a surface differing only far above its minimum.
        grid = np.linspace(-np.pi, np.pi, 200)
        shape = 30.0 * (1.0 - np.cos(grid))          # 0 at the well, 60 at the top
        moved = shape.copy()
        moved[shape > DRIFT_CEILING_KJMOL] += 5.0    # only the barrier moves

        difference = np.abs(shape - moved)
        judged = shape <= DRIFT_CEILING_KJMOL

        assert difference.max() > SETTLED_DRIFT_KJMOL      # over the whole grid
        assert difference[judged].max() <= SETTLED_DRIFT_KJMOL  # where it counts


class TestARefusedRunStillHandsBackWhatItHas:
    """The refusal says "the surface below is a snapshot of that filling" and
    then returned None, so a run whose wells had converged to half of RT left
    nothing to plot. "Not converged" and "not available" are different
    claims."""

    def _refused(self, tmp_path):
        grid = _grid()
        hills = _hills_reaching(_double_well(grid), grid)
        return compute_surface(_written(tmp_path, hills), np.tile([0., 1.], 20))

    def test_the_snapshot_comes_back(self, tmp_path) -> None:
        result = self._refused(tmp_path)
        assert result["refused"]
        assert result["surface"] is not None
        assert len(result["surface"]) == len(result["grid"])

    def test_it_is_marked_provisional(self, tmp_path) -> None:
        """Both are worth having; only one is worth quoting."""
        result = self._refused(tmp_path)
        assert result["provisional"] is True

    def test_the_flag_distinguishes_the_two(self, tmp_path) -> None:
        """`provisional` is what tells a reader whether the numbers beside it
        are an answer or a progress report, so it has to be present either
        way rather than only when it is true."""
        result = self._refused(tmp_path)
        assert "provisional" in result
        assert isinstance(result["provisional"], bool)

        import inspect

        from fastmdxplora.simulation import metad_surface

        source = inspect.getsource(metad_surface.compute_surface)
        # Set on the accepted path too, not only the refused one.
        assert source.count('"provisional"') >= 2

    def test_a_coordinate_that_never_moved_still_has_nothing(
        self, tmp_path
    ) -> None:
        """The one refusal where None is the honest answer: no surface exists
        along a coordinate that did not move."""
        still = Hills(np.arange(50.0), np.full(50, 0.3), np.full(50, SIGMA),
                      np.full(50, 1.0), GAMMA)
        result = compute_surface(_written(tmp_path, still), np.full(50, 0.3))
        assert result["surface"] is None
