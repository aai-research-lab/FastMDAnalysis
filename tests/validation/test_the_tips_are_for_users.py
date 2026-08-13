"""Every tip the GUI shows is written for the person choosing, not for us.

The analysis picker surfaces each module docstring's first two paragraphs.
Three modules opened instead with the maintainer's story of the bug that
motivated them -- "written to `pmf.json` and stopped there: no figure, no
entry in the analysis manifest" -- and the GUI rendered that history,
faithfully, beside a checkbox where somebody was deciding whether to run a
PMF. Accurate, and useless to them.

The history is worth keeping; it now lives below the two-paragraph window
the GUI reads. What this test enforces is the window itself: everything a
user is shown -- analysis summaries and details, analysis option help,
every schema field's help -- answers a user's questions (what is this, how
do I read it, when does it apply) rather than narrating the software's.

The markers are the shared fingerprints of the narrative voice, chosen
tight: "used to" and "was written" appear in perfectly good user guidance
("the seed is written into the manifest"), so only phrases that cannot be
anything but history are refused.
"""

from __future__ import annotations

import pytest

from fastmdxplora.analysis.describe import describe_all, explain_analysis
from fastmdxplora.config.schema import PHASE_SCHEMAS, TOP_LEVEL

#: Phrases that only ever narrate the software's past, never guide a user.
NARRATIVE = (
    "stopped there",
    "no mention in the report",
    "existed to produce",
    "each produced a curve",
    "the same gap",
    "until now",
    "this session",
    "cost somebody",
    "no manifest entry",
)


def _narrative_in(text: str) -> list[str]:
    lowered = (text or "").lower()
    return [marker for marker in NARRATIVE if marker in lowered]


class TestAnalysisTips:
    @pytest.mark.parametrize("name", sorted(describe_all()))
    def test_the_tip_answers_a_user(self, name: str) -> None:
        doc = explain_analysis(name)
        assert doc["summary"], f"{name} has no summary to show"
        assert doc["detail"], (
            f"{name} has no detail paragraph: the GUI would show its title "
            "twice and explain nothing")
        assert len(doc["detail"]) > 80, (
            f"{name}'s detail is too thin to help anyone choose")
        for part in ("summary", "detail"):
            found = _narrative_in(doc[part])
            assert not found, (
                f"{name}'s {part} narrates the software's history "
                f"({found}) where a user needs to know what the analysis "
                "does. Move the story below the docstring's second "
                "paragraph; the GUI reads only the first two.")

    @pytest.mark.parametrize("name", sorted(describe_all()))
    def test_option_help_guides_rather_than_narrates(self, name: str) -> None:
        for option in describe_all()[name]:
            found = _narrative_in(option.help or "")
            assert not found, (
                f"{name}.{option.name}'s help narrates history: {found}")


class TestSettingTips:
    def test_every_field_help_is_for_users(self) -> None:
        offenders = []
        for phase, group in list(PHASE_SCHEMAS.items()) + [("top", TOP_LEVEL)]:
            for field in group.fields:
                found = _narrative_in(field.help or "")
                if found:
                    offenders.append(f"{phase}.{field.name}: {found}")
        assert not offenders, offenders


class TestTheThreeThatTaughtThis:
    """The original offenders, pinned by what they must now say."""

    def test_pmf_tells_a_user_how_to_read_a_profile(self) -> None:
        detail = explain_analysis("pmf")["detail"]
        assert "Minima" in detail and "kJ/mol" in detail
        assert "overlap" in detail  # the honesty a user must know about

    def test_metad_surface_explains_the_landscape_and_provisional(self) -> None:
        detail = explain_analysis("metad_surface")["detail"]
        assert "basins" in detail.lower()
        assert "provisional" in detail

    def test_steered_work_says_what_it_is_not(self) -> None:
        summary_and_detail = " ".join(explain_analysis("steered_work").values())
        assert "smooth" in summary_and_detail
        # The free-energy caveat sits in paragraph three, past the GUI's
        # window on purpose -- but the curve-reading guidance must be in it.
        assert "snapped past" in summary_and_detail

    def test_the_history_is_kept_below_the_window(self) -> None:
        """Demoted, not deleted: the docstrings still carry their story for
        maintainers, past the two paragraphs the GUI reads."""
        import fastmdxplora.analysis.pmf as pmf_module

        assert "stopped there" in (pmf_module.__doc__ or "")
