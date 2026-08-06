"""Ten analyses averaged over the whole production run without asking whether
it had settled, or how many independent samples the average rested on.

Both follow from one quantity, so both are fixed by one piece of arithmetic:
the statistical inefficiency is the number of frames per independent sample,
which says where the relaxation ended and what the mean is worth.
"""

from __future__ import annotations

import numpy as np
import pytest

from fastmdxplora.analysis.equilibration import (
    MINIMUM_EFFECTIVE_SAMPLES,
    correlation_is_resolved,
    detect_equilibration,
    statistical_inefficiency,
    summarise,
)


def _correlated(phi, n=200000, seed=0):
    """An AR(1) series, whose statistical inefficiency is known exactly.

    ``x[t] = phi x[t-1] + noise`` has autocorrelation ``phi**t``, so
    ``g = 1 + 2 sum phi**t = (1 + phi) / (1 - phi)``. That makes it the one
    case where the answer can be checked rather than eyeballed.
    """
    rng = np.random.RandomState(seed)
    values = np.zeros(n)
    for i in range(1, n):
        values[i] = phi * values[i - 1] + rng.normal()
    return values


class TestCountingIndependentSamples:
    """A trajectory written every picosecond, from a system decorrelating over
    a hundred, has a hundred times fewer independent samples than frames."""

    @pytest.mark.parametrize("phi", [0.0, 0.5, 0.8, 0.9])
    def test_it_matches_the_known_answer(self, phi) -> None:
        expected = (1.0 + phi) / (1.0 - phi)
        measured = statistical_inefficiency(_correlated(phi))
        assert measured == pytest.approx(expected, rel=0.05)

    def test_uncorrelated_frames_each_count_once(self) -> None:
        rng = np.random.RandomState(0)
        assert statistical_inefficiency(rng.normal(size=50000)) == \
            pytest.approx(1.0, abs=0.05)

    def test_a_constant_series_has_nothing_to_correlate(self) -> None:
        """No fluctuations, so no correlation time describes them. One rather
        than a division by zero."""
        assert statistical_inefficiency(np.full(500, 2.5)) == 1.0

    def test_too_few_frames_to_ask(self) -> None:
        assert statistical_inefficiency(np.array([1.0, 2.0])) == 1.0

    def test_it_never_reports_fewer_frames_than_samples(self) -> None:
        """An inefficiency below one would mean frames carrying more than
        their own information."""
        rng = np.random.RandomState(3)
        for series in (rng.normal(size=1000), _correlated(0.7, 5000, seed=4)):
            assert statistical_inefficiency(series) >= 1.0


class TestFindingWhereARunSettled:
    def test_a_relaxation_is_discarded(self) -> None:
        """Several relaxation times of it, which is what leaves the average
        describing the equilibrium rather than the approach to it."""
        rng = np.random.RandomState(0)
        tau = 200.0
        frames = np.arange(6000)
        series = 3.0 * np.exp(-frames / tau) + rng.normal(0, 0.15, frames.size)

        discard, _g, _n = detect_equilibration(series)

        assert 2 * tau < discard < 6 * tau

    def test_an_already_settled_run_keeps_its_frames(self) -> None:
        rng = np.random.RandomState(0)
        discard, _g, effective = detect_equilibration(rng.normal(size=5000))

        assert discard < 500
        assert effective > 4000

    def test_the_answer_never_rests_on_the_last_third(self) -> None:
        """Whatever the series does, an average of its tail is not an average
        of the run."""
        rng = np.random.RandomState(0)
        awkward = np.concatenate([rng.normal(size=3000),
                                  rng.normal(5.0, 0.01, 2000)])
        discard, _g, _n = detect_equilibration(awkward)
        assert discard <= int(awkward.size * 2 / 3)

    def test_a_short_series_is_left_alone(self) -> None:
        assert detect_equilibration(np.arange(5.0))[0] == 0


class TestAnErrorBarThatMeansSomething:
    def test_correlation_widens_it(self) -> None:
        """Computed as though every frame were independent, the error on a
        strongly correlated series is understated several-fold -- which is how
        a difference between two systems becomes significant on paper without
        being real."""
        series = _correlated(0.95, 20000, seed=1)
        settled, refusal = summarise(series)

        assert refusal is None
        naive = float(np.std(series, ddof=1) / np.sqrt(series.size))
        assert settled.standard_error > 4 * naive

    def test_the_spread_is_not_the_error(self) -> None:
        """One is a property of the system, the other of how long it was
        watched, and reporting either for the other is a common way to make a
        result look tighter or looser than it is."""
        settled, _ = summarise(_correlated(0.8, 20000, seed=2))
        assert settled.standard_deviation > settled.standard_error

    def test_a_run_with_too_few_independent_samples_is_refused(self) -> None:
        """Uncorrelated, so the count is measurable -- and there are still not
        enough of them."""
        rng = np.random.RandomState(0)
        settled, refusal = summarise(rng.normal(size=8))

        assert refusal is not None
        assert "independent samples" in refusal
        assert settled is not None, "the numbers are still returned to be read"

    def test_the_refusal_says_recording_more_often_will_not_help(self) -> None:
        """It is the commonest wrong response: the frames are correlated, so
        more of them are more of the same. Said where the count is measurable
        and small, since a run that cannot measure its own correlation gets
        the more specific refusal instead."""
        rng = np.random.RandomState(0)
        _settled, refusal = summarise(rng.normal(size=8))
        assert "recording them more often will not help" in refusal

    def test_length_is_not_the_test(self) -> None:
        """A long correlated run can hold fewer independent samples than a
        short uncorrelated one."""
        rng = np.random.RandomState(6)
        long_and_correlated = _correlated(0.999, 4000, seed=7)
        short_and_clean = rng.normal(size=200)

        _s1, refused = summarise(long_and_correlated)
        _s2, accepted = summarise(short_and_clean)

        assert refused is not None and accepted is None
        assert "not long against its own correlation time" in refused


class TestARunTooShortToMeasureItsOwnCorrelation:
    """The estimate of the correlation is itself unreliable when the run is
    short against the thing being estimated -- and it fails in the flattering
    direction, reporting more independent samples than there are.

    An AR(1) series with a true inefficiency of 2000 gave 361 from four
    thousand frames. The independent-sample count built from it said eleven,
    cleared a threshold of ten, and the truth was two.
    """

    def test_a_series_that_cannot_resolve_its_correlation_is_caught(
        self
    ) -> None:
        assert not correlation_is_resolved(_correlated(0.999, 4000, seed=7))

    def test_one_that_can_is_not(self) -> None:
        assert correlation_is_resolved(_correlated(0.99, 20000, seed=7))
        assert correlation_is_resolved(_correlated(0.9, 20000, seed=7))

    def test_even_a_long_run_can_fail_it(self) -> None:
        """Twenty thousand frames still cannot resolve a correlation time of
        a thousand: the estimate comes back 1129 against a true 1999, and the
        check says so rather than reporting the sample count as measured."""
        assert not correlation_is_resolved(_correlated(0.999, 20000, seed=7))

    def test_asking_where_the_correlation_decayed_does_not_work(self) -> None:
        """Worth pinning, because it is the obvious check. A sample
        autocorrelation sums to about -1/2 whatever the series, so it goes
        negative on its own: on the unresolvable series it crosses zero at a
        lag where the real correlation is still 0.7, and a check built on
        that crossing passes the run it is meant to catch.
        """
        series = _correlated(0.999, 4000, seed=7)
        fluctuation = series - series.mean()
        variance = float(np.mean(fluctuation ** 2))

        crossing = next(
            lag for lag in range(1, series.size - 1)
            if np.mean(fluctuation[:series.size - lag] * fluctuation[lag:])
            / variance <= 0.0)
        true_correlation_there = 0.999 ** crossing

        assert crossing < 500
        assert true_correlation_there > 0.5, (
            "the correlation is still strong where the estimator says it ended")

    def test_the_refusal_says_the_count_is_an_upper_bound(self) -> None:
        _settled, refusal = summarise(_correlated(0.999, 4000, seed=7))
        assert "upper bound" in refusal
        assert "longer run is the only remedy" in refusal

    def test_a_short_series_with_a_correlation_cannot_resolve_it(self) -> None:
        assert not correlation_is_resolved(np.arange(30.0))

    def test_but_frames_near_independence_have_nothing_to_resolve(self) -> None:
        """A short uncorrelated series has no correlation time for a longer
        run to pin down, and putting one through the halving comparison only
        measures the noise in the comparison: discarding a single frame
        flipped the verdict on twenty-five."""
        rng = np.random.RandomState(0)
        assert correlation_is_resolved(rng.normal(size=25))
        assert correlation_is_resolved(rng.normal(size=24))

    def test_nothing_to_average_is_said_plainly(self) -> None:
        settled, refusal = summarise(np.array([1.0]))
        assert settled is None
        assert "nothing to average" in refusal

    def test_frames_that_are_not_numbers_are_dropped_first(self) -> None:
        rng = np.random.RandomState(0)
        series = rng.normal(size=5000)
        series[[10, 20, 30]] = np.nan

        settled, refusal = summarise(series)
        assert refusal is None
        assert np.isfinite(settled.mean) and np.isfinite(settled.standard_error)

    def test_the_record_carries_every_number_behind_the_mean(self) -> None:
        """So a reader can check the claim rather than take it."""
        settled, _ = summarise(_correlated(0.5, 20000, seed=8))
        assert set(settled.as_record()) == {
            "discard", "statistical_inefficiency", "effective_samples",
            "mean", "standard_error", "standard_deviation"}

    def test_the_threshold_is_a_judgement_a_study_can_make(self) -> None:
        rng = np.random.RandomState(0)
        marginal = rng.normal(size=40)

        _s, strict = summarise(marginal, minimum_effective_samples=1000)
        _s, lenient = summarise(marginal, minimum_effective_samples=5)

        assert strict is not None and lenient is None
        assert MINIMUM_EFFECTIVE_SAMPLES == 10.0
