"""AUD40: the growth loop used the cube's arithmetic on every box.

`_solvate_with_room_for_the_cutoff` exists to notice a box narrower than
twice the cutoff and re-solvate with more padding. It computed the extra
padding as `shortfall / 2`, which is the relation for a cube: a cube's
narrowest periodic width is `maxSize + 2 * padding`, so 0.30 nm of
shortfall wants 0.15 nm more padding.

A dodecahedron's narrowest width is that quantity over sqrt(2). The same
padding buys sqrt(2) less width, and the loop undershot every non-cubic
box by exactly that factor.

Measured on the run it killed. V3, alanine dipeptide, 1.20 nm padding, 1.0
nm cutoff, dodecahedron: the box came out 1.70 nm across, the loop grew to
1.45 nm intending 2.20, and produced 2.05 against a floor of 2.00. NPT
then contracted the box 1.9% toward real water's density -- which is the
barostat's whole job -- and the run died on the second move with "The
periodic box size has decreased to less than twice the nonbonded cutoff."

Two things were wrong and the second is why the first mattered: the
increment used the wrong factor, and the target was the bare floor rather
than the floor plus the contraction the next phase was about to apply.
"""

from __future__ import annotations

import math

import pytest

from fastmdxplora.setup.prepare import (
    NARROWEST_WIDTH_PER_SIZE,
    NPT_CONTRACTION_MARGIN,
    padding_that_reaches,
)


def _width_after_growing(smallest: float, padding: float, grown: float,
                         shape: str) -> float:
    """d(width)/d(padding) = 2f, which is the relation under test."""
    return smallest + 2.0 * NARROWEST_WIDTH_PER_SIZE[shape] * (grown - padding)


class TestTheFactorsAreTheBoxVectors:
    """Not chosen: read off the matrix OpenMM builds and constrains."""

    def test_they_match_openmm_s_own_vectors(self) -> None:
        pytest.importorskip("openmm")

        expected = {
            "cube": min(1.0, 1.0, 1.0),
            "dodecahedron": min(1.0, 1.0, math.sqrt(2) / 2),
            "octahedron": min(1.0, 2 * math.sqrt(2) / 3, math.sqrt(6) / 3),
        }
        for shape, want in expected.items():
            assert NARROWEST_WIDTH_PER_SIZE[shape] == pytest.approx(want)

    def test_a_dodecahedron_is_the_narrow_one(self) -> None:
        assert (NARROWEST_WIDTH_PER_SIZE["dodecahedron"]
                < NARROWEST_WIDTH_PER_SIZE["octahedron"]
                < NARROWEST_WIDTH_PER_SIZE["cube"])


class TestEveryShapeReachesItsTarget:

    @pytest.mark.parametrize("shape", sorted(NARROWEST_WIDTH_PER_SIZE))
    def test_the_grown_box_clears_the_floor_with_the_margin(
            self, shape: str) -> None:
        cutoff, smallest, padding = 1.0, 1.70, 1.20
        grown = padding_that_reaches(
            smallest_nm=smallest, padding_nm=padding,
            nonbonded_cutoff_nm=cutoff, box_shape=shape)
        width = _width_after_growing(smallest, padding, grown, shape)
        assert width == pytest.approx(2.0 * cutoff * NPT_CONTRACTION_MARGIN)

    def test_a_box_that_already_clears_is_not_grown(self) -> None:
        assert padding_that_reaches(
            smallest_nm=3.0, padding_nm=1.2, nonbonded_cutoff_nm=1.0,
            box_shape="dodecahedron") == 1.2

    def test_an_unknown_shape_is_treated_as_a_cube(self) -> None:
        """Conservative in the wrong direction is still wrong, but a shape
        this does not know about should not silently get a factor invented
        for it."""
        assert padding_that_reaches(
            smallest_nm=1.7, padding_nm=1.2, nonbonded_cutoff_nm=1.0,
            box_shape="something-else") == padding_that_reaches(
            smallest_nm=1.7, padding_nm=1.2, nonbonded_cutoff_nm=1.0,
            box_shape="cube")


class TestTheRunThisKilled:
    """V3, and the numbers are the ones the run logged."""

    LOGGED_SMALLEST = 1.70
    LOGGED_PADDING = 1.20
    CUTOFF = 1.0
    CONTRACTION = 0.98085   # 0.9408 -> 0.997 g/mL, linear

    def test_the_old_arithmetic_undershot_by_root_two(self) -> None:
        old_grown = self.LOGGED_PADDING + (2.0 * self.CUTOFF
                                           - self.LOGGED_SMALLEST) / 2.0 + 0.1
        old_width = _width_after_growing(
            self.LOGGED_SMALLEST, self.LOGGED_PADDING, old_grown,
            "dodecahedron")
        # 2.05 against the 2.0527 the run actually built.
        assert old_width == pytest.approx(2.05, abs=0.01)
        assert old_width < 2.0 * self.CUTOFF * NPT_CONTRACTION_MARGIN

    def test_the_old_box_could_not_survive_the_contraction(self) -> None:
        """This is why it is a defect rather than a tight fit: no
        equilibration of that box to correct density was possible."""
        old_grown = self.LOGGED_PADDING + (2.0 * self.CUTOFF
                                           - self.LOGGED_SMALLEST) / 2.0 + 0.1
        old_width = _width_after_growing(
            self.LOGGED_SMALLEST, self.LOGGED_PADDING, old_grown,
            "dodecahedron")
        assert old_width * self.CONTRACTION < 2.0 * self.CUTOFF + 0.02

    def test_the_new_box_survives_it(self) -> None:
        grown = padding_that_reaches(
            smallest_nm=self.LOGGED_SMALLEST, padding_nm=self.LOGGED_PADDING,
            nonbonded_cutoff_nm=self.CUTOFF, box_shape="dodecahedron")
        width = _width_after_growing(
            self.LOGGED_SMALLEST, self.LOGGED_PADDING, grown, "dodecahedron")
        assert width * self.CONTRACTION > 2.0 * self.CUTOFF

    def test_the_margin_covers_the_measured_contraction(self) -> None:
        """2-3.5% linear on the systems measured. A margin that did not
        cover it would leave the same defect with a bigger number."""
        assert NPT_CONTRACTION_MARGIN > 1.035
