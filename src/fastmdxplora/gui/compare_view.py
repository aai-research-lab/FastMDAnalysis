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
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

__all__ = ["compared"]


def compared(first: Path | str, second: Path | str) -> dict[str, Any]:
    """Everything the Compare page draws for two studies or runs."""
    from fastmdxplora.gui.workspace import studies_compared

    a, b = Path(first), Path(second)
    said = studies_compared(a, b)
    if not said.get("ok"):
        return said
    said["first"] = {**said["first"], **_where(a)}
    said["second"] = {**said["second"], **_where(b)}
    said["series"] = _series_of_both(a, b)
    return said


def _where(folder: Path) -> dict[str, Any]:
    """A study's name as the page says it: a run by its study and by what
    sets it apart among that study's runs ("chignolin, seed 1")."""
    from fastmdxplora.gui.exploration import study_a_run_belongs_to

    try:
        run = study_a_run_belongs_to(folder)
    except Exception:  # noqa: BLE001 - named by its folder instead
        run = None
    if not run:
        return {"said": folder.name, "run_of": None}
    study = Path(run["path"]).name
    return {"said": f"{study}, {run['label']}", "run_of": {"path": run["path"], "name": study,
                                                          "label": run["label"]}}


def _series_of_both(a: Path, b: Path) -> list[dict[str, Any]]:
    """Each series over time that both recorded, in the order the Analysis
    page lists them, with the numbers of each."""
    from fastmdxplora.gui.series import series_over_time, series_payload

    try:
        mine = series_over_time(a)["series"]
        theirs = {entry["analysis"] for entry in series_over_time(b)["series"]}
    except Exception:  # noqa: BLE001 - the page stands without them
        return []
    both: list[dict[str, Any]] = []
    for entry in mine:
        name = entry["analysis"]
        if name not in theirs:
            continue
        one, other = series_payload(a, name), series_payload(b, name)
        if not (one.get("ok") and other.get("ok")
                and one.get("kind") == "time" and other.get("kind") == "time"):
            continue
        both.append({"analysis": name, "label": entry["label"], "unit": entry["unit"],
                     "first": _drawn(one), "second": _drawn(other)})
    return both


def _drawn(series: dict[str, Any]) -> dict[str, Any]:
    """What a line on the page needs of a series: its points, what its x
    runs over, and the mean with where its averaging starts."""
    return {"x": series["x"], "y": series["y"], "x_label": series["x_label"],
            "thinned": series.get("thinned"), "mean": series.get("mean")}
