"""The comparison script, now that it lives in the repository.

It cannot be run here: it needs ProLIF and MDAnalysis, which the package
does not depend on, and a finished run to read. What can be checked is that
it still imports, that the thresholds it was pre-registered with are the
ones in the file, and that the lessons it was hardened with have not been
tidied out of it -- each of those was a cluster-side failure that cost an
afternoon, and a rewrite that loses one will cost the next afternoon too.
"""

from __future__ import annotations

import inspect

import pytest


@pytest.fixture(scope="module")
def module():
    return pytest.importorskip(
        "fastmdxplora.validation.cross_tool",
        reason="the cross-tool comparison needs ProLIF and MDAnalysis")


class TestThePreRegisteredThresholds:
    """Changed after the fact, a pre-registration is not one.

    These are the numbers the protocol was written with before any
    comparison had been run, so they are checked rather than trusted.
    """

    def test_the_three_tolerances_are_unchanged(self, module):
        assert module.HARMONIZED_OCCUPANCY_TOL_PP == 5.0
        assert module.OBSERVABLE_TOL_NM == 1e-3
        assert module.NEGATIVE_OCCUPANCY_MAX_PCT == 5.0


class TestTheLessonsAreStillInIt:
    """Each of these was found by a failure on the cluster."""

    def test_it_still_heals_the_periodic_boundary(self, module):
        source = inspect.getsource(module)
        assert "image_molecules" in source, (
            "a reference tool reading raw coordinates sees molecules split "
            "across the box, and distance-guessed bonds cannot heal them")

    def test_it_still_re_roots_cluster_paths(self, module):
        source = inspect.getsource(module)
        assert "fallback" in source.lower()

    def test_the_ceiling_groups_by_protein_atom(self, module):
        """Summing correlated pairs convicted a passing negative control.

        Ring carbons touching one protein atom fire in the same frames, so
        their occupancies add to something the residue never had.
        """
        source = inspect.getsource(module)
        assert "by_patom" in source or "protein_atom" in source

    def test_the_negative_control_is_judged_on_the_floor(self, module):
        source = inspect.getsource(module)
        assert "FLOOR" in source or "floor" in source
