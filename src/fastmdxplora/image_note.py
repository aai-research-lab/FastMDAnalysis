"""How a mean's image note is said, read from its record: one reader for the
report, every page that shows the mean, and the Agent.

The image check (:mod:`fastmdxplora.analysis.image_check`) finds how close
the solute came to its own periodic images in each frame, and every mean of
those frames carries a note of it where the solute came within the run's
cutoff of its image in the frames the mean uses (``solute_image`` on the
mean's record, its numbers; the analysis's ``findings.solute_image``, its
paragraph). The note does not change the mean's verdict: a mean determined
is still determined, its error kept, and the note is said beside it.

The report said the note and the Overview, the Analysis page's table and the
Agent did not, so one mean read "may be biased" on one page and "Determined"
alone on the next. Every reader now reads the note here, from the record:
:func:`note_of` (the mean's numbers), :func:`said_short` (beside the
verdict, "may be biased: 54 of 270 frames used came within 1 nm of the
solute's own image"), :func:`paragraph_of` (the analysis's paragraph, one
click or one look away), :func:`finding_said` (the image check's own
finding, in a sentence, for its row and its card) and :func:`longer_run_said`
(beside the advice of a longer run, for a run in contact with its image in
every frame).

Light on purpose: no analysis code is imported, so a page reads it without
loading the trajectory tools.
"""

from __future__ import annotations

from typing import Any

#: What a note's count is of, as it says it: the trajectory's frames a mean
#: uses, the energy rows a thermodynamic mean averages (each read at the
#: trajectory frame nearest its time), or, for an analysis without a main
#: mean, the frames its means use.
COUNTED_OVER = {
    "frames": "frames used",
    "energy rows": "energy rows used",
    "means' frames": "frames its means use",
}


def share_said(part: int, whole: int) -> str:
    """A share of frames as a percentage, never rounded to 0% or 100%."""
    value = 100.0 * part / whole if whole else 0.0
    if 0 < value < 1:
        return f"{value:.2g}%"
    if 99 < value < 100:
        return f"{value:.1f}%"
    return f"{value:.0f}%"


def _count(value: Any) -> int | None:
    return value if isinstance(value, int) and not isinstance(value, bool) else None


def note_of(record: Any) -> dict[str, Any] | None:
    """The note a mean record (or an analysis's findings, or a reweighted
    average) carries: its ``solute_image`` where the solute came within the
    cutoff of its image in some of the frames it uses; else None. The image
    check's own finding is not a note (it has no ``frames_used``)."""
    note = record.get("solute_image") if isinstance(record, dict) else None
    if not isinstance(note, dict):
        return None
    used, under = _count(note.get("frames_used")), _count(note.get("frames_under_cutoff"))
    if not (used and under and used > 0 and under > 0):
        return None
    return note


def said_short(note: dict[str, Any] | None) -> str:
    """The note beside a mean's verdict, in a line: "may be biased: 54 of
    270 frames used (20%) came within 1 nm of the solute's own image". An
    empty string where there is no note."""
    if not note_of({"solute_image": note}):
        return ""
    assert note is not None
    used, under = int(note["frames_used"]), int(note["frames_under_cutoff"])
    cutoff = note.get("cutoff_nm")
    within = (f"within {float(cutoff):g} nm of" if isinstance(cutoff, (int, float))
              and not isinstance(cutoff, bool) else "within the cutoff of")
    over = str(note.get("counted_over") or "frames")
    noun = COUNTED_OVER.get(over, "frames used")
    share = share_said(under, used)
    if over == "energy rows":
        return (f"may be biased: in {under:,} of {used:,} {noun} ({share}) the solute "
                f"was {within} its own image")
    return (f"may be biased: {under:,} of {used:,} {noun} ({share}) came {within} the "
            "solute's own image")


def paragraph_of(findings: Any) -> str:
    """The analysis's own paragraph of its note (``findings.solute_image``,
    counted over the frames its main mean uses, or its means'), or of the
    image check's finding; an empty string where there is none."""
    image = findings.get("solute_image") if isinstance(findings, dict) else None
    if not isinstance(image, dict):
        return ""
    return str(image.get("said") or "")


def longer_run_said(finding: Any) -> str:
    """Where the image check found the solute in contact with its own image
    (under 0.5 nm) in every frame it read, that a longer run, which would
    determine the means, does not lift that bias, and what does; an empty
    string otherwise. One sentence for the report, the Overview's lead and
    the Agent's answers, which advised only a longer run."""
    if not isinstance(finding, dict) or "frames_used" in finding:
        return ""
    frames, close = _count(finding.get("frames")), _count(finding.get("frames_under_close"))
    if not frames or close != frames:
        return ""
    return ("A longer run would determine the means; it would not lift the bias of the "
            "solute's own image, in contact with it in every frame: a longer run in the "
            "same box samples the same contact, and a larger box, in a new study, is what "
            "lifts it.")


def recorded_finding(root: Any) -> dict[str, Any] | None:
    """The image check's finding as a run recorded it (its
    ``analysis/image_check/options.json``, read through `record_access`);
    None where there is none, or it cannot be read."""
    import json
    from pathlib import Path

    from fastmdxplora import record_access

    try:
        record = json.loads(record_access.read_text(
            Path(root) / "analysis" / "image_check" / "options.json", encoding="utf-8"))
    except (OSError, ValueError, TypeError):
        return None
    found = (record.get("findings") or {}).get("solute_image") if isinstance(record, dict) \
        else None
    return found if isinstance(found, dict) else None


def finding_said(finding: Any) -> str:
    """The image check's finding in a sentence, its first: how close the
    solute came to its own image, and in how many frames it was within the
    cutoff and in contact. For the check's own row and card; the rest of the
    paragraph (the closest frame and atoms, the cutoff's source, the solute
    and how it was grouped) is one click away."""
    if not isinstance(finding, dict) or "frames_used" in finding:
        return ""
    said = " ".join(str(finding.get("said") or "").split())
    cut = said.find(". ")
    return said[:cut + 1] if cut > 0 else said
