"""g(r) reported a hydration shell where the curve was flat.

The tallest point of a curve always exists. On a 0.1 ns Trp-cage run the
curve spanned 0.000 to 1.074 and its tallest point sat at 1.782 nm -- the
last bin before half the box -- which was recorded as
`first_peak_nm: 1.7818`, `first_peak_height: 1.074`, with no qualification.
A reader plots that and reports a hydration shell at 1.78 nm.

The same machinery gives the T4 radius of gyration a `not_a_measurement`
verdict naming its own unreliability. This brings g(r) into line: two
criteria, both read off the curve rather than set as thresholds.
"""
import numpy as np
import pytest

from fastmdxplora.analysis.rdf import _is_that_a_peak


RADII = np.linspace(0.003, 1.794, 358)


def _peak_of(g):
    return int(np.argmax(np.where(RADII > 0.15, g, -np.inf)))


def _curve(peak_nm=None, height=0.0, noise=0.02):
    """A curve with a chosen peak and a chosen amount of scatter.

    The scatter is a deterministic sawtooth rather than a random draw. With
    random noise these assertions turn on where one realization happens to
    fall relative to the criterion -- an early draft passed on one seed and
    failed on another at the same nominal noise level, which is a flaky test
    dressed as a physics test.
    """
    wobble = noise * np.sign(np.sin(np.arange(RADII.size) * 2.0))
    g = 1.0 + wobble
    if peak_nm is not None:
        g = g + height * np.exp(-((RADII - peak_nm) / 0.03) ** 2)
    return np.where(RADII < 0.15, 0.0, g)


def test_a_genuine_hydration_shell_is_reported():
    g = _curve(peak_nm=0.28, height=1.4)
    assert _is_that_a_peak(RADII, g, _peak_of(g)) == ""


def test_a_weak_peak_is_still_a_peak_when_the_curve_is_clean():
    """Height alone does not decide it. A small peak on a quiet curve is
    resolved; the same peak on a noisy one is not."""
    g = _curve(peak_nm=0.30, height=0.30, noise=0.01)
    assert _is_that_a_peak(RADII, g, _peak_of(g)) == ""


def test_the_same_peak_buried_in_noise_is_refused():
    g = _curve(peak_nm=0.30, height=0.30, noise=0.25)
    verdict = _is_that_a_peak(RADII, g, _peak_of(g))
    assert "not distinguishable from the noise" in verdict


def test_a_flat_curve_reports_no_peak():
    g = _curve(noise=0.025)
    verdict = _is_that_a_peak(RADII, g, _peak_of(g))
    assert verdict
    assert "longer run" in verdict or "edge of what" in verdict


def test_the_tallest_point_at_the_edge_is_the_box_not_a_shell():
    """A first peak is a local maximum. The last bins are where the
    periodic cell runs out."""
    g = _curve(noise=0.005)
    g[-1] = 3.0                      # tall, but at the boundary
    verdict = _is_that_a_peak(RADII, g, _peak_of(g))
    assert "half the box" in verdict


def test_the_verdict_says_what_would_change_it():
    """A refusal that does not say what to do is half a refusal."""
    for g in (_curve(noise=0.025), _curve(peak_nm=0.3, height=0.3, noise=0.25)):
        verdict = _is_that_a_peak(RADII, g, _peak_of(g))
        assert any(word in verdict for word in ("longer run", "larger box"))


def test_a_curve_too_short_to_judge_is_not_second_guessed():
    """With almost no bulk region there is no noise estimate, and inventing
    one would be worse than reporting the number."""
    short_r = np.linspace(0.2, 0.4, 5)
    short_g = np.array([1.0, 1.2, 2.0, 1.1, 1.0])
    assert _is_that_a_peak(short_r, short_g, 2) == ""
