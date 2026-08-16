"""The harness that measures the guardrails, and its own honesty."""

from __future__ import annotations

import pytest

from fastmdxplora.validation.corpus import (
    CLEAN,
    DEFECTS,
    Case,
    run_case,
    run_corpus,
)


class TestTheCorpusItself:
    def test_every_defect_is_detected(self):
        result = run_corpus(DEFECTS, CLEAN)
        assert result["defects"]["detection_rate"] == 1.0, (
            f"missed: {[m['case'] for m in result['defects']['missed']]}")

    def test_nothing_fires_on_ordinary_work(self):
        result = run_corpus(DEFECTS, CLEAN)
        assert result["clean"]["false_refusal_rate"] == 0.0, (
            f"refused: {[f['case'] for f in result['clean']['false_refusals']]}")
        assert not result["clean"]["other_disagreements"]

    def test_both_corpora_are_populated(self):
        """A detection rate over an empty clean corpus is not a measurement.

        The rate that makes the other one mean anything is the false-refusal
        rate, and it needs cases where the right answer is silence.
        """
        assert len(DEFECTS) >= 5
        assert len(CLEAN) >= 3

    def test_every_case_records_why(self):
        for case in DEFECTS + CLEAN:
            assert case.because, f"{case.name} does not say why"
            assert case.expect in ("refused", "qualified", "proceeded")


class TestTheHarnessCanFail:
    """A measurement that cannot come out badly is not measuring.

    These give the harness cases whose answers are wrong on purpose, so a
    passing corpus above means the software behaved rather than the
    harness being unable to say otherwise.
    """

    def test_a_guardrail_that_does_not_fire_is_recorded_as_missed(self):
        silent = Case("a defect nothing catches", lambda: {"value": 1.0},
                      "refused", "it should have been refused")
        outcome = run_case(silent)
        assert outcome.observed == "proceeded"
        assert not outcome.agreed

    def test_refusing_ordinary_work_is_recorded_as_a_false_refusal(self):
        def _refuses():
            raise ValueError("no")

        result = run_corpus(
            [], [Case("ordinary work", _refuses, "proceeded", "it is fine")])
        assert result["clean"]["false_refusal_rate"] == 1.0

    def test_refusing_for_the_wrong_reason_does_not_count(self):
        """Otherwise a guardrail could pass by failing at something else.

        A case that crashes on a typo in its own fixture would read as a
        detection, and the corpus would report a rate it had not earned.
        """
        def _wrong_reason():
            raise ValueError("the file was not found")

        case = Case("a defect", _wrong_reason, "refused",
                    "it should be refused for the stated reason",
                    mentioning="periodic boundary")
        outcome = run_case(case)
        assert outcome.observed == "refused"
        assert not outcome.agreed
        assert "not for the stated reason" in outcome.detail

    def test_a_qualification_is_neither_a_refusal_nor_a_clean_pass(self):
        """Collapsing the three would hide what the software is built on."""
        qualified = Case(
            "answered with a caveat",
            lambda: {"rdf": {"capped": "stopped at half the box"}},
            "qualified", "the number stands with a statement attached")
        assert run_case(qualified).agreed

        assert run_case(Case(
            "answered with a caveat", qualified.run, "proceeded",
            "wrongly expected to be clean")).agreed is False
        assert run_case(Case(
            "answered with a caveat", qualified.run, "refused",
            "wrongly expected to be refused")).agreed is False

    def test_a_qualification_one_level_in_is_still_seen(self):
        """An analysis writes its findings under its own name.

        Reading only the top level scored a guardrail that had fired as a
        miss, which is the instrument failing rather than the software.
        """
        nested = Case(
            "nested", lambda: {"thermo": {"not_a_measurement": "fixed box"}},
            "qualified", "the caveat is one level in")
        assert run_case(nested).agreed
