"""Two studies, or two runs, compared on one page.

Comparison was in three places, each with its own UI: two studies chosen in
All studies drew their settings and means under the cards, a study of
several runs drew its runs' table and overlay on its Analysis page, and the
Viewer played another study beside the one open. This is what the one
Compare page reads: the settings in which the two differ and the means each
recorded, a difference marked only where it is resolved (as
:func:`fastmdxplora.gui.workspace.studies_compared` says it), each series
both recorded over time to overlay, and where each is, so a run is named by
its study and by what sets it apart. The Viewer plays the second beside the
first from it.

Each of the two is read held to itself (`server._confined`, through
`record_access`): a record linked out of it, or to a file of it that is not
one of its results, is not read, and the page says so, never taking it for
a record the study does not have.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Callable

__all__ = ["compared"]


def compared(first: Path | str, second: Path | str) -> dict[str, Any]:
    """Everything the Compare page draws for two studies or runs."""
    from fastmdxplora.gui.workspace import sides_compared

    a, b = Path(first), Path(second)
    sides = [_side(a), _side(b)]
    for side, folder in zip(sides, (a, b)):
        if not side.get("ok"):
            return {"ok": False, "reason": _refused_why(side, folder)}
    said = {"ok": True, **sides_compared(sides[0]["compared"], sides[1]["compared"])}
    places = [_where(a), _where(b)]
    said["first"] = {**sides[0]["head"], **places[0]["where"]}
    said["second"] = {**sides[1]["head"], **places[1]["where"]}
    series, drawn = _series_of_both(a, b, sides)
    said["series"] = series
    not_read = [line for line in (
        _not_read_said(said["first"]["said"], sides[0], drawn[0]),
        _not_read_said(said["second"]["said"], sides[1], drawn[1]),
        _not_read_said(places[0]["study"], places[0]),
        _not_read_said(places[1]["study"], places[1]),
    ) if line]
    if not_read:
        said["not_read"] = " ".join(not_read)
    return said


def _held(folder: Path, read: Callable[[Path], dict[str, Any]]) -> dict[str, Any]:
    """``read(folder)`` held to the folder's own records (`server._confined`)."""
    from fastmdxplora.gui.server import _confined

    return _confined(folder, read)


def _side(folder: Path) -> dict[str, Any]:
    """What the comparison reads of one of the two, held to it: what it is,
    its Config, the phases and analyses it ran, its means, and the series it
    recorded over time."""
    from fastmdxplora.gui.series import series_over_time
    from fastmdxplora.gui.workspace import compared_side, head_of
    from fastmdxplora.refusals import refusal_of

    def read(inside: Path) -> dict[str, Any]:
        try:
            side = compared_side(inside)
        except Exception as exc:  # noqa: BLE001 - said, not raised
            return {"ok": False, "reason": refusal_of(exc).message}
        try:
            series = series_over_time(inside)["series"]
        except Exception:  # noqa: BLE001 - the page stands without them
            series = []
        # Named as it was opened, not as its links resolve; the folder they
        # resolve to too, as a study opened is said (`/api/explore/switch`).
        return {"ok": True, "compared": side,
                "head": {**head_of(inside), "path": str(folder), "resolved": str(inside)},
                "series": series}

    return _held(folder, read)


def _refused_why(side: dict[str, Any], folder: Path) -> str:
    """Why one of the two could not be compared: a record of it not read is
    said so, named with the study, never as a study without it."""
    said = _not_read_said(folder.name, side)
    return said or str(side.get("reason") or "They could not be compared.")


def _not_read_said(name: str, *answers: dict[str, Any]) -> str:
    """The records the held reads of one study refused, each with why, in a
    sentence beginning with the study's name: "" where every one was read."""
    from fastmdxplora.gui.server import _refusals_said

    names: list[str] = []
    why: dict[str, str] = {}
    for answer in answers:
        for refused in answer.get("unavailable") or ():
            if refused not in names:
                names.append(refused)
        why.update(answer.get("unavailable_why") or {})
    said = _refusals_said({"unavailable": sorted(names), "unavailable_why": why})
    return f"{name}: {said}" if said else ""


def _where(folder: Path) -> dict[str, Any]:
    """A study's name as the page says it: a run by its study and by what
    sets it apart among that study's runs ("chignolin, seed 1"), read from
    its study's records held to that study."""
    from fastmdxplora.gui.exploration import study_a_run_belongs_to

    alone = {"where": {"said": folder.name, "run_of": None}, "study": folder.name}
    if folder.parent.name != "runs":
        return alone
    study = folder.parent.parent

    def read(inside: Path) -> dict[str, Any]:
        try:
            return {"run": study_a_run_belongs_to(inside / "runs" / folder.name)}
        except Exception:  # noqa: BLE001 - named by its folder instead
            return {"run": None}

    answer = _held(study, read)
    run = answer.get("run")
    if not run:
        return {**answer, **alone, "study": study.name}
    said = f"{study.name}, {run['label']}"
    return {**answer, "study": study.name,
            "where": {"said": said, "run_of": {"path": str(study), "name": study.name,
                                               "label": run["label"]}}}


def _series_of_both(a: Path, b: Path,
                    sides: list[dict[str, Any]]) -> tuple[list[dict[str, Any]], list[dict]]:
    """Each series over time that both recorded, in the order the Analysis
    page lists them, with the numbers of each; and the held reads of each
    of the two that drew them."""
    from fastmdxplora.gui.series import series_payload

    theirs = {entry["analysis"] for entry in sides[1]["series"]}
    common = [entry for entry in sides[0]["series"] if entry["analysis"] in theirs]
    if not common:
        return [], [{}, {}]

    def payloads(inside: Path) -> dict[str, Any]:
        return {"payloads": {entry["analysis"]: series_payload(inside, entry["analysis"])
                             for entry in common}}

    drawn = [_held(a, payloads), _held(b, payloads)]
    both: list[dict[str, Any]] = []
    for entry in common:
        name = entry["analysis"]
        one, other = (side.get("payloads", {}).get(name) or {} for side in drawn)
        if not (one.get("ok") and other.get("ok")
                and one.get("kind") == "time" and other.get("kind") == "time"):
            continue
        both.append({"analysis": name, "label": entry["label"], "unit": entry["unit"],
                     "first": _drawn(one), "second": _drawn(other)})
    return both, drawn


def _drawn(series: dict[str, Any]) -> dict[str, Any]:
    """What a line on the page needs of a series: its points, what its x
    runs over, and the mean with where its averaging starts."""
    return {"x": series["x"], "y": series["y"], "x_label": series["x_label"],
            "thinned": series.get("thinned"), "mean": series.get("mean")}
