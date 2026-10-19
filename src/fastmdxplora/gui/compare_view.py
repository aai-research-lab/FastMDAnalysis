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

import os
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
    places = [_where(a), _where(b)]
    replicas = _as_replicas(a, b, places, sides)
    said = {"ok": True, **sides_compared(sides[0]["compared"], sides[1]["compared"],
                                         replicas.get("errors"))}
    if replicas.get("errors"):
        _the_replicas_said(said, replicas)
    elif replicas.get("unfinished"):
        names = " and ".join(places[i]["where"]["said"] for i in replicas["unfinished"])
        said["replicas_not_used"] = (
            f"These are runs that differ only by random seed, but {names} did not finish: the "
            "spread of the finished runs' means is the error of a run that went its full "
            "length, not of one that ended early. Each difference is set against the errors "
            "the two runs estimated for themselves, which replicas can show too tight.")
    elif replicas.get("unknown"):
        said["replicas_unknown"] = (
            f"Whether these two runs are replicas is not known: {places[0]['study']}'s "
            "batch_manifest.json could not be read. Each difference is set against the errors "
            "the two runs estimated for themselves, which replicas can show too tight.")
    said["first"] = {**sides[0]["head"], **places[0]["where"]}
    said["second"] = {**sides[1]["head"], **places[1]["where"]}
    series, drawn, unread = _series_of_both(a, b, sides)
    said["series"] = series
    if unread:
        said["series_not_read"] = unread
    # Two runs of one study: its records' sentence said once.
    studies = [(places[0]["study_path"], places[0]), (places[1]["study_path"], places[1])]
    if studies[0][0] and studies[0][0] == studies[1][0]:
        studies = studies[:1]
    not_read = [line for line in (
        _not_read_said(said["first"]["said"], sides[0], drawn[0]),
        _not_read_said(said["second"]["said"], sides[1], drawn[1]),
        *(_not_read_said(place["study"], place, replicas if index == 0 else {})
          for index, (_, place) in enumerate(studies)),
    ) if line]
    outside = replicas.get("outside") or []
    if outside:
        not_read.append(f"{places[0]['study']}: " + ", ".join(outside) + ": "
                        + ("its" if len(outside) == 1 else "their")
                        + " manifest places " + ("it" if len(outside) == 1 else "them")
                        + " outside the study, so " + ("it was" if len(outside) == 1 else "they were")
                        + " not read and not counted among the replicas.")
    if not_read:
        said["not_read"] = " ".join(not_read)
    return said


def _as_replicas(a: Path, b: Path, places: list[dict[str, Any]],
                 sides: list[dict[str, Any]]) -> dict[str, Any]:
    """Two runs of one study whose runs are replicas (they differ only by
    seed, `batch.aggregate`): the error of each one's mean, by quantity,
    as the replicas give it (`batch.aggregate.replica_error`: the spread of
    the means of the study's finished runs, or the run's own error where
    that is larger), read held to the study. Two repeats of one
    measurement differ by its error, which the spread of the repeats
    measures and each run's own estimate can understate.

    With ``errors``: ``runs``, the finished runs; ``counts``, by quantity,
    how many of them recorded the mean its spread is taken over; and
    ``kinds``, by quantity, which error each of the two has ("spread" or
    "own", where the run's own is larger). ``unfinished``: the two are
    replicas but one of them did not finish, so the spread is not used.
    ``outside``: finished runs its manifest places outside the study, left
    out unread (named by their place, never looked for). Empty for two
    studies, for two runs of a sweep of anything else, and where the
    study's records were not read."""
    study = places[0].get("study_path")
    if not study or study != places[1].get("study_path"):
        return {}
    wanted = [name for name, row in _both_measured(sides).items() if row]
    if not wanted:
        return {}

    def read(inside: Path) -> dict[str, Any]:
        from fastmdxplora.batch.aggregate import finished_members, replica_error
        from fastmdxplora.gui.workspace import means_of

        try:
            members, replicas, why = finished_members(inside)
        except (OSError, ValueError):
            # Whether they are repeats of one measurement is not known.
            return {"unknown": True}
        if not replicas:
            return {}
        # Each finished run by its place in the study, without following a
        # link or looking outside it: one its manifest places elsewhere is
        # named, and left out.
        kept: list[Path] = []
        outside: list[str] = []
        for run_id, member in members:
            try:
                kept.append(Path(os.path.abspath(member)).relative_to(inside))
            except ValueError:
                outside.append(f"runs/{run_id or member.name}")
        finished = {run_id for run_id, _ in members}
        unfinished = [i for i, run in enumerate((a, b)) if run.name not in finished]
        if unfinished:
            # A run that stopped or failed is a shorter measurement than
            # those the spread is taken over: their spread is not its error.
            return {"unfinished": unfinished, "outside": outside}
        means = [means_of(inside / place) for place in kept]
        errors: dict[str, tuple[float, float]] = {}
        kinds: dict[str, tuple[str, str]] = {}
        counts: dict[str, int] = {}
        for name in wanted:
            values = [(m.get(name) or {}).get("mean") for m in means]
            values = [v for v in values if v is not None]
            own = [sides[i]["compared"]["means"][name]["error"] for i in (0, 1)]
            pair = [replica_error(values, own[i]) for i in (0, 1)]
            spread = replica_error(values, None)
            if None not in pair and spread is not None:
                errors[name] = (pair[0], pair[1])
                kinds[name] = tuple("own" if own[i] is not None and own[i] > spread else "spread"
                                    for i in (0, 1))
                counts[name] = len(values)
        return {"errors": errors, "kinds": kinds, "counts": counts, "runs": len(members),
                "why": why, "outside": outside}

    from fastmdxplora.gui.server import _confined

    return _confined(Path(study), read)


def _the_replicas_said(said: dict[str, Any], replicas: dict[str, Any]) -> None:
    """Each difference between two replicas marked with the error it was
    set against: the replicas' spread where that is each one's error, each
    run's own where that was the larger, and how many of the runs' means
    the spread was taken over."""
    kinds, counts = replicas.get("kinds") or {}, replicas.get("counts") or {}
    for row in said.get("measures") or ():
        versus, name = row.get("versus"), row.get("analysis")
        if not versus or name not in kinds:
            continue
        mine = list(kinds[name])
        versus["error_kinds"] = mine
        versus["error_of"] = ("replicas" if mine == ["spread", "spread"]
                              else "own" if mine == ["own", "own"] else "mixed")
        versus["from_runs"] = counts.get(name)
    said["replicas"] = {"runs": replicas["runs"], "why": replicas["why"]}


def _both_measured(sides: list[dict[str, Any]]) -> dict[str, bool]:
    """Each quantity, and whether both of the two have its mean with an
    error, the only means compared."""
    def usable(side: dict[str, Any] | None) -> bool:
        return bool(side) and side.get("mean") is not None and side.get("error") is not None

    mine, theirs = (side["compared"]["means"] for side in sides)
    return {name: usable(mine.get(name)) and usable(theirs.get(name))
            for name in set(mine) | set(theirs)}


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

    alone = {"where": {"said": folder.name, "run_of": None}, "study": folder.name,
             "study_path": None}
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
        return {**answer, **alone, "study": study.name, "study_path": str(study)}
    said = f"{study.name}, {run['label']}"
    return {**answer, "study": study.name, "study_path": str(study),
            "where": {"said": said, "run_of": {"path": str(study), "name": study.name,
                                               "label": run["label"]}}}


def _series_of_both(a: Path, b: Path, sides: list[dict[str, Any]]
                    ) -> tuple[list[dict[str, Any]], list[dict], list[str]]:
    """Each series over time that both recorded, in the order the Analysis
    page lists them, then production's energies; the held reads of each of
    the two that drew them; and the series both recorded that were not
    drawn because a record of one was not read."""
    from fastmdxplora.gui.series import series_payload

    theirs = {entry["analysis"] for entry in sides[1]["series"]}
    common = [entry for entry in sides[0]["series"] if entry["analysis"] in theirs]

    def payloads(inside: Path) -> dict[str, Any]:
        return {"payloads": {entry["analysis"]: series_payload(inside, entry["analysis"])
                             for entry in common},
                "energies": _energies(inside)}

    drawn = [_held(a, payloads), _held(b, payloads)]
    both: list[dict[str, Any]] = []
    unread: list[str] = []
    for entry in common:
        name = entry["analysis"]
        one, other = (side.get("payloads", {}).get(name) or {} for side in drawn)
        if one.get("unavailable") or other.get("unavailable"):
            unread.append(entry["label"])
            continue
        if not (one.get("ok") and other.get("ok")
                and one.get("kind") == "time" and other.get("kind") == "time"):
            continue
        both.append({"analysis": name, "label": entry["label"], "unit": entry["unit"],
                     "first": _drawn(one), "second": _drawn(other)})
    from fastmdxplora.gui.overview_view import THERMODYNAMICS

    mine, others = (side.get("energies") or {} for side in drawn)
    for key, label, unit in THERMODYNAMICS:
        if key in mine and key in others:
            both.append({"analysis": key, "label": label, "unit": unit,
                         "first": mine[key], "second": others[key]})
    return both, drawn, unread


def _energies(root: Path) -> dict[str, dict[str, Any]]:
    """Production's thermodynamic quantities over its time, from OpenMM's
    energy file as the run's Overview reads it for their means (its pieces
    joined where it was carried on); none where it was not read."""
    from fastmdxplora.gui.overview_view import THERMODYNAMICS, NotRead, _production_energies
    from fastmdxplora.gui.series import MOST_POINTS

    try:
        rows = _production_energies(root)
    except Exception:  # noqa: BLE001 - the page stands without them
        return {}
    if not rows or isinstance(rows, NotRead):
        return {}
    found: dict[str, dict[str, Any]] = {}
    for key, _, _ in THERMODYNAMICS:
        points = [(t, v) for t, v in ((_number(row.get("simulation_time_ns")),
                                       _number(row.get(key))) for row in rows)
                  if t is not None and v is not None]
        if len(points) < 2:
            continue
        every = max(1, -(-len(points) // MOST_POINTS))
        kept = points[::every]
        # Production's clock, from zero where production began (the runner
        # sets it so), named as every other series on the card names it
        # (`series.analysed_axis`), so the page says them in one unit.
        found[key] = {"x": [t for t, _ in kept], "y": [v for _, v in kept],
                      "x_label": "Time (ns)",
                      "thinned": every if every > 1 else None, "mean": None}
    return found


def _number(value: Any) -> float | None:
    import math

    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _drawn(series: dict[str, Any]) -> dict[str, Any]:
    """What a line on the page needs of a series: its points, what its x
    runs over, and the mean with where its averaging starts."""
    return {"x": series["x"], "y": series["y"], "x_label": series["x_label"],
            "thinned": series.get("thinned"), "mean": series.get("mean")}
