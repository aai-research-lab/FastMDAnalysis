"""What the Viewer's timeline marks on the frames of the run it plays.

The timeline under the Viewer's molecule (viewer-timeline.js) marks the
moments a person may want to go to: each cluster's medoid, the frame its
clustering points to as the one most like the rest of its cluster; and each
time a study run until it is determined judged it (`stopping.json`), at the
production it had reached. The series under it mark where each analysis
began averaging; those come with the series (`/api/series`).

Only records are read: a medoid's frame as the clustering wrote it in its
populations table, a judgement's production as the stopping record keeps
it. Frames are the trajectory's own numbers, as `/api/frames-info` gives
them; a time is production's clock.
"""

from __future__ import annotations

import csv
import math
import re
from pathlib import Path
from typing import Any

__all__ = ["timeline_payload"]

#: The clustering's tables of populations, one for each method it ran.
_POPULATIONS = re.compile(r"cluster_(?P<method>[a-z0-9_]+)_populations\.csv")


def _number(value: Any) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _medoids(root: Path) -> list[dict[str, Any]]:
    marks: list[dict[str, Any]] = []
    folder = root / "analysis" / "cluster"
    if not folder.is_dir():
        return marks
    for table in sorted(folder.glob("cluster_*_populations.csv")):
        named = _POPULATIONS.fullmatch(table.name)
        if not named:
            continue
        method = named.group("method").replace("_", " ")
        try:
            with table.open(encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle))
        except (OSError, csv.Error):
            continue
        for row in rows:
            frame = _number(row.get("medoid_frame"))
            share = _number(row.get("fraction"))
            if frame is None:
                continue
            marks.append({
                "kind": "medoid",
                "frame": int(frame),
                "time_ns": _number(row.get("medoid_time_ns")),
                "label": (f"The medoid of {method} cluster {row.get('cluster', '?')}"
                          + (f", {share:.0%} of the frames" if share is not None else "")),
            })
    return marks


def _judgements(root: Path) -> list[dict[str, Any]]:
    from fastmdxplora.gui.stopping_view import _record_of

    found = _record_of(root)
    if not found:
        return []
    rounds = found["record"].get("rounds") or []
    marks = []
    for number, entry in enumerate(rounds, start=1):
        if not isinstance(entry, dict):
            continue
        at = _number(entry.get("production_ns"))
        if at is None:
            continue
        decision = str(entry.get("decision") or entry.get("outcome") or "").replace("_", " ")
        marks.append({"kind": "check", "frame": None, "time_ns": at,
                      "label": f"Judged after round {number}" + (f": {decision}" if decision else "")})
    return marks


def timeline_payload(root: str | Path) -> dict[str, Any]:
    """The marks of the run at ``root``, each with its frame or its time."""
    base = Path(root)
    return {"ok": True, "marks": _medoids(base) + _judgements(base)}
