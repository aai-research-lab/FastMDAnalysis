"""Checks on the corpus itself, which need no network.

A corpus that claims coverage it does not have is worse than none: it
reads as evidence. These run on every ordinary `pytest`.
"""

from __future__ import annotations

import pytest

from .structures import CORPUS, kinds


class TestTheCorpusIsHonest:
    @pytest.mark.parametrize("expectation", CORPUS, ids=lambda e: e.label)
    def test_a_refusal_or_caveat_says_what_is_wrong(self, expectation) -> None:
        """A refusal and a caveat are different: one is a structure
        turned away, the other a build that is correct as far as it
        goes. Conflating them made the corpus assert that 5WYZ fails,
        which it does not."""
        for recorded in (expectation.refused, expectation.caveat):
            if recorded:
                # An explanation rather than a label.
                assert len(recorded) > 80, expectation.label
        assert not (expectation.refused and expectation.caveat), (
            f"{expectation.label} is both refused and caveated")

    def test_the_kinds_are_distinct(self) -> None:
        assert len(set(kinds())) == len(CORPUS)

    def test_the_corpus_covers_the_classes_it_claims_to(self) -> None:
        covered = {expectation.kind for expectation in CORPUS}
        for wanted in (
            "soluble", "soluble-ligand", "multimer-glycoprotein",
            "cofactor-metal", "membrane-helical", "membrane-barrel",
            "membrane-fusion", "membrane-large", "nucleic-acid",
            "protein-nucleic",
        ):
            assert wanted in covered, wanted
