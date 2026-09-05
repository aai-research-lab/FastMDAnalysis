"""The harness that measures the guardrails, and its own honesty."""

from __future__ import annotations


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


class TestTheCleanCorpusIsLargeEnoughToClaimSomething:
    """Zero false refusals in seven cases excludes very little.

    The exact 95% upper bound on a rate given zero failures in n trials is
    1 - 0.05**(1/n): 35% at seven, 18% at fifteen. The registration asked
    for fifteen because specificity is the half of V6 that is not a
    foregone conclusion -- a checker that refused everything would score
    perfectly on the defect corpus, and only the clean corpus separates
    the two.
    """

    def test_the_clean_corpus_meets_its_registration(self) -> None:
        from fastmdxplora.validation.corpus import CLEAN

        assert len(CLEAN) >= 15, (
            f"{len(CLEAN)} clean cases; zero refusals there is consistent "
            f"with a true rate below {100 * (1 - 0.05 ** (1 / len(CLEAN))):.0f}%, "
            "and the pre-registration named fifteen")

    def test_every_case_is_named_once(self) -> None:
        """Two cases sharing a name make the rates uninterpretable."""
        from fastmdxplora.validation.corpus import CLEAN, DEFECTS

        names = [c.name for c in DEFECTS] + [c.name for c in CLEAN]
        assert len(names) == len(set(names))

    def test_no_case_is_written_and_left_unwired(self) -> None:
        """`_a_torsion_read_as_a_straight_line` was defined, complete, and
        referenced nowhere -- a guardrail case that exists and never runs
        is the same as one that does not exist. Wiring it in showed it does
        not fire on its own data, so it stays out and stays named here
        rather than sitting silently in the module."""
        import inspect

        from fastmdxplora.validation import corpus

        wired = {c.run for c in corpus.DEFECTS} | {c.run for c in corpus.CLEAN}
        known_unwired = {"_a_torsion_read_as_a_straight_line"}
        defined = {
            name for name, obj in vars(corpus).items()
            if name.startswith("_") and inspect.isfunction(obj)
            and obj.__module__ == corpus.__name__
            and not name.startswith("__")
        }
        helpers = {"_classify", "_numpy", "_hills", "_peptide", "_gas",
                   "_qualification"}
        orphans = {
            n for n in defined - helpers
            if vars(corpus)[n] not in wired
            and inspect.signature(vars(corpus)[n]).parameters == {}
        } - known_unwired
        assert not orphans, (
            f"case functions defined and never run: {sorted(orphans)}")
