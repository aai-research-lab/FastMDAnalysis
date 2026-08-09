"""The averages of a biased run, said in the order they should be read.

The analyses of a metadynamics trajectory each report a mean, and every one
of those means is an average over a distribution the bias flattened on
purpose. The reweighted value is the measurement of the system; the raw one
is a property of the sampling. So the reweighted value leads here and the raw
one follows it, rather than the arrangement the rest of the report would
naturally produce, where the corrected number appears in a section of its own
and the uncorrected one sits under every analysis heading above it.

The effective sample size is printed beside the number and not in a footnote.
A reweighted mean over five hundred frames whose weight sits in seven of them
is a mean over seven, and there is no arrangement of a document in which that
should be readable without the seven.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

#: Below this the correction is reported with the caution stated in the same
#: breath as the number, rather than only in the notes underneath.
THIN_EFFECTIVE_FRAMES = 50.0


def _record_path(project_root: Path) -> Path:
    return (Path(project_root) / "analysis" / "reweighted"
            / "reweighted_averages.json")


def load_reweighted(project_root: Path) -> dict[str, Any] | None:
    """What the analysis phase computed, or None where no bias applied.

    Absent for plain MD, for a steered pull and for an umbrella window, and
    that absence is not a gap: none of those has a set of weights that
    recovers an equilibrium average, which the methods text already says.
    """
    import json

    path = _record_path(project_root)
    if not path.is_file():
        return None
    try:
        record = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    if not isinstance(record, dict) or not record.get("quantities"):
        return None
    return record


def _number(value: Any, digits: int = 4) -> str:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return "—"
    if number != number:  # nan
        return "—"
    return f"{number:.{digits}g}"


def _with_spread(mean: Any, spread: Any) -> str:
    text = _number(mean)
    spread_text = _number(spread, 2)
    if text == "—" or spread_text == "—":
        return text
    return f"{text} ± {spread_text}"


def reweighted_line(record: dict[str, Any] | None, analysis: str) -> str | None:
    """One line for an analysis's own section, leading with the correction.

    A reader who goes straight to the RMSD heading should not have to find
    another section to learn that the number above it is an average over a
    biased ensemble.
    """
    if not record:
        return None
    for item in record.get("quantities") or []:
        if item.get("analysis") != analysis:
            continue
        ess = record.get("effective_sample_size")
        frames = record.get("n_frames")
        return (
            f"**Reweighted mean: {_with_spread(item.get('reweighted_mean'), item.get('reweighted_std'))}**"
            f" — the equilibrium average, recovered from the deposited bias on "
            f"{_number(ess, 3)} effective frames of {frames}. "
            f"The average over the biased frames themselves is "
            f"{_number(item.get('raw_mean'))}"
            f" ({_signed(item.get('shift_percent'))})."
        )
    return None


def _percent(value: Any) -> str:
    """Signed, because the direction the bias moved a number is the point."""
    try:
        number = float(value)
    except (TypeError, ValueError):
        return "—"
    if number != number:
        return "—"
    return f"{number:+.1f}%"


def _signed(value: Any) -> str:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return "change not available"
    if number != number:
        return "change not available"
    return f"{number:+.1f}% from reweighting"


def reweighted_section(project_root: Path,
                       record: dict[str, Any] | None = None) -> str:
    """The corrected averages as a table, with what they rest on."""
    if record is None:
        record = load_reweighted(project_root)
    if not record:
        return ""

    ess = record.get("effective_sample_size")
    frames = record.get("n_frames")
    settled = bool(record.get("settled"))

    lines = ["### Averages after reweighting", ""]
    lines.append(
        "Production was biased, so the trajectory is not a Boltzmann "
        "ensemble and an average over its frames is not a measurement of the "
        "system. Each frame is weighted by exp((V − c(t))/RT), where V is the "
        "bias it was actually sampled under and c(t) is the Tiwary–Parrinello "
        "offset, which recovers the unbiased average. **The reweighted column "
        "is the result; the biased column is shown so the size of the "
        "correction is visible.**")
    lines.append("")

    caution = ""
    if isinstance(ess, (int, float)) and ess < THIN_EFFECTIVE_FRAMES:
        caution = (
            " These averages therefore rest on very few independent frames, "
            "and the spreads beside them are correspondingly wide.")
    lines.append(
        f"The weights carry {_number(ess, 3)} effective frames of {frames}."
        + caution)
    if not settled:
        lines.append("")
        lines.append(
            "_The bias had not settled when the run ended, so these are "
            "approximate._")
    lines.append("")

    lines.append("| Quantity | Reweighted (equilibrium) | Biased trajectory | Change |")
    lines.append("| --- | --- | --- | --- |")
    for item in record.get("quantities") or []:
        lines.append(
            f"| {item.get('label', item.get('analysis', '—'))} "
            f"| {_with_spread(item.get('reweighted_mean'), item.get('reweighted_std'))} "
            f"| {_with_spread(item.get('raw_mean'), item.get('raw_std'))} "
            f"| {_percent(item.get('shift_percent'))} |")
    lines.append("")

    figure = (Path(project_root) / "analysis" / "reweighted"
              / "reweighted_averages.png")
    if figure.is_file():
        lines.append("![Effect of reweighting]"
                     "(analysis/reweighted/reweighted_averages.png)")
        lines.append("")

    warnings = record.get("warnings") or []
    if warnings:
        for warning in warnings:
            lines.append(f"- _{warning}_")
        lines.append("")

    # Said plainly rather than left to be inferred from the table's absence.
    lines.append(
        "_Analyses whose result is not one value per frame — the clustering "
        "and the dimensionality reduction — are not reweighted, and the "
        "populations and projections they report are those of the biased "
        "ensemble._")
    lines.append("")
    return "\n".join(lines)
