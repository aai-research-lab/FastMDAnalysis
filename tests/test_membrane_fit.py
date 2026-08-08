"""Fitting the bilayer, rather than assuming one and checking afterwards.

Orienting by principal axes asks which way a structure is longest, which is
the wrong question wherever a soluble domain carries the shape: 6B73's
longest axis runs through its G protein. And every check of a supplied
orientation compares populations across the whole molecule, so a T4L fusion
or a beta-barrel's polar lumen swamps the belt being looked for -- which is
why the calibration set overlapped, with a porin refused at 1.17 and
adenylate kinase refused at 0.91.

These tests exercise the arithmetic, which needs neither OpenMM nor a real
structure. The measurement against known proteins is a separate exercise.
"""

from __future__ import annotations

import numpy as np
import pytest

from fastmdxplora.setup.membrane_fit import (
    fit_membrane_slab,
    hemisphere_directions,
    rotation_onto_z,
)


def _belted(normal, count=1200, radius=1.6, half=1.5, seed=0):
    """A cylinder with an apolar band around its middle and polar ends."""
    rng = np.random.default_rng(seed)
    normal = np.asarray(normal, dtype=float)
    normal = normal / np.linalg.norm(normal)
    across = np.cross(normal, [1.0, 0.0, 0.0])
    if np.linalg.norm(across) < 1e-6:
        across = np.cross(normal, [0.0, 1.0, 0.0])
    across = across / np.linalg.norm(across)
    other = np.cross(normal, across)

    along = rng.uniform(-2.0 * half, 2.0 * half, count)
    angle = rng.uniform(0.0, 2.0 * np.pi, count)
    points = (np.outer(along, normal)
              + radius * (np.outer(np.cos(angle), across)
                          + np.outer(np.sin(angle), other)))
    inside = np.abs(along) < half
    return points, np.where(inside, 0.5, 0.02), np.where(inside, 0.02, 0.5)


def _globular(count=1200, radius=2.2, apolar_fraction=0.35, seed=1):
    """A soluble protein: exposed surface mixed, no band anywhere."""
    rng = np.random.default_rng(seed)
    direction = rng.normal(size=(count, 3))
    direction /= np.linalg.norm(direction, axis=1)[:, None]
    points = direction * radius
    apolar = rng.random(count) < apolar_fraction
    return points, np.where(apolar, 0.5, 0.0), np.where(apolar, 0.0, 0.5)


class TestTheSearchCoversDirections:
    def test_they_are_unit_vectors(self) -> None:
        directions = hemisphere_directions(500)
        assert np.allclose(np.linalg.norm(directions, axis=1), 1.0)

    def test_a_hemisphere_not_a_sphere(self) -> None:
        """A slab is unchanged by flipping its normal, so searching both
        halves does the same work twice."""
        assert (hemisphere_directions(500)[:, 2] >= 0).all()

    def test_they_are_spread_rather_than_bunched(self) -> None:
        directions = hemisphere_directions(500)
        separations = directions @ directions.T
        np.fill_diagonal(separations, -1.0)
        # No two directions closer than a few degrees.
        assert separations.max() < 0.999


class TestTheFitFindsTheMembrane:
    @pytest.mark.parametrize(
        "normal", [(0, 0, 1), (1, 0, 0), (0.3, -0.5, 0.8), (-0.6, 0.6, 0.5)])
    def test_the_normal_is_recovered_from_any_frame(self, normal) -> None:
        """The fit decides which way is up; it does not need to be told."""
        fit = fit_membrane_slab(*_belted(normal))
        wanted = np.asarray(normal, dtype=float)
        wanted /= np.linalg.norm(wanted)
        assert abs(float(np.dot(fit.normal, wanted))) > 0.98

    def test_the_thickness_is_the_one_that_was_there(self) -> None:
        fit = fit_membrane_slab(*_belted((0, 0, 1), half=1.5))
        assert fit.thickness_nm == pytest.approx(3.0, abs=0.3)

    def test_the_slab_is_mostly_apolar(self) -> None:
        fit = fit_membrane_slab(*_belted((0, 0, 1)))
        assert fit.hydrophobic_fraction > 0.9

    def test_a_thinner_belt_is_found_thinner(self) -> None:
        thin = fit_membrane_slab(*_belted((0, 0, 1), half=1.1))
        thick = fit_membrane_slab(*_belted((0, 0, 1), half=1.9))
        assert thin.thickness_nm < thick.thickness_nm


class TestASolubleProteinDoesNotFit:
    def test_the_best_slab_still_scores_badly(self) -> None:
        """There is always a best slab. What distinguishes a membrane protein
        is that its best slab buries more apolar surface than polar, and a
        soluble protein's does not."""
        fit = fit_membrane_slab(*_globular())
        assert fit.score_nm2 < 0

    def test_and_covers_an_ordinary_surface(self) -> None:
        fit = fit_membrane_slab(*_globular(apolar_fraction=0.35))
        assert fit.hydrophobic_fraction < 0.55

    def test_the_two_are_far_apart(self) -> None:
        """The measure they replace put a porin and a kinase within 0.3 of
        each other, on the wrong sides of the line."""
        membrane = fit_membrane_slab(*_belted((0, 0, 1)))
        soluble = fit_membrane_slab(*_globular())
        assert membrane.score_nm2 > 0 > soluble.score_nm2
        assert membrane.hydrophobic_fraction - soluble.hydrophobic_fraction > 0.4


class TestTheFitRefusesToGuess:
    def test_too_few_atoms(self) -> None:
        points, apolar, polar = _belted((0, 0, 1), count=10)
        assert fit_membrane_slab(points, apolar, polar) is None

    def test_mismatched_arrays(self) -> None:
        points, apolar, polar = _belted((0, 0, 1))
        assert fit_membrane_slab(points, apolar[:-5], polar) is None

    def test_no_surface_at_all(self) -> None:
        points, apolar, polar = _belted((0, 0, 1))
        assert fit_membrane_slab(points, apolar * 0, polar * 0) is None


class TestTheRotation:
    @pytest.mark.parametrize(
        "normal", [(0, 0, 1), (0, 0, -1), (1, 0, 0), (0.3, -0.5, 0.8)])
    def test_it_puts_the_normal_on_z(self, normal) -> None:
        rotation = rotation_onto_z(normal)
        wanted = np.asarray(normal, dtype=float)
        wanted /= np.linalg.norm(wanted)
        assert np.allclose(rotation @ wanted, [0.0, 0.0, 1.0], atol=1e-9)

    def test_it_is_a_rotation(self) -> None:
        """Not a reflection: a mirrored protein is a different molecule."""
        for normal in [(0, 0, 1), (1, 0, 0), (0.3, -0.5, 0.8), (-0.2, 0.9, 0.1)]:
            rotation = rotation_onto_z(normal)
            assert np.allclose(rotation @ rotation.T, np.eye(3), atol=1e-9)
            assert float(np.linalg.det(rotation)) == pytest.approx(1.0)

    def test_the_fit_and_the_rotation_agree(self) -> None:
        points, apolar, polar = _belted((0.3, -0.5, 0.8))
        fit = fit_membrane_slab(points, apolar, polar)
        rotated = points @ rotation_onto_z(fit.normal).T
        # After rotating, the apolar band lies about the fitted centre in z.
        band = rotated[apolar > 0.1][:, 2]
        assert abs(float(np.median(band)) - fit.centre_nm) < 0.4
