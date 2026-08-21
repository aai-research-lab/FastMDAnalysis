"""An analysis that refuses says why, where the person is looking.

`rdf` on a default configuration cannot succeed: its default `selection_b`
is `water and name O` and the simulation phase's default `save_selection`
excludes water. It raised with a message that names the problem exactly --
"selection_b 'water and name O' matched no atoms, so there is no
distribution between the two groups to report" -- and that message went to
the debug log while the console showed `✗ error` and nothing else. The
remedy, `save_selection: all`, was printed by a different phase.

The message was always carried on AnalysisResult.message. The row renderer
had no parameter for it and the caller dropped it.
"""
import io
from contextlib import redirect_stdout

import pytest

from fastmdxplora.utils.presenter import SessionPresenter


def row(**kwargs) -> str:
    buf = io.StringIO()
    with redirect_stdout(buf):
        SessionPresenter().analysis_table_row(**kwargs)
    return buf.getvalue()


RDF_MESSAGE = (
    "rdf: selection_b 'water and name O' matched no atoms, so there is no "
    "distribution between the two groups to report."
)


def test_the_reason_is_printed_not_only_logged():
    out = row(name="rdf", status="error", path="analysis/rdf/", elapsed=0.0,
              name_width=16, reason=RDF_MESSAGE)
    assert "matched no atoms" in out
    assert "selection_b" in out


def test_the_analysis_name_is_not_repeated():
    """It is already the first column."""
    out = row(name="rdf", status="error", path="analysis/rdf/", elapsed=0.0,
              name_width=16, reason=RDF_MESSAGE)
    # once in the row itself, not again at the head of the reason
    assert out.count("rdf") == 2, out          # name column + path
    assert "rdf: selection_b" not in out


def test_a_successful_analysis_gains_nothing():
    out = row(name="rmsd", status="ok", path="analysis/rmsd/", elapsed=0.34,
              name_width=16)
    assert out.count("\n") == 1


@pytest.mark.parametrize("useless", ["KeyError", "IndexError", "", "   ", None])
def test_a_message_that_explains_nothing_is_left_out(useless):
    """A bare exception class occupies the line where an explanation
    belongs while saying only that something raised. The traceback is in
    the log for whoever needs it."""
    out = row(name="cluster", status="error", path="analysis/cluster/",
              elapsed=0.0, name_width=16,
              reason=f"cluster: {useless}" if useless else useless)
    assert out.count("\n") == 1, out


def test_a_long_reason_wraps_under_the_row():
    out = row(name="rdf", status="error", path="analysis/rdf/", elapsed=0.0,
              name_width=16, reason=RDF_MESSAGE)
    lines = out.rstrip("\n").split("\n")
    assert len(lines) > 1
    assert all(len(line) < 200 for line in lines)
