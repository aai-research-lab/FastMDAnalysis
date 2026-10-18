"""Tags and a note a person gives a study, kept in the study's folder.

A study's records say what it is and what it found; nothing said what it
means to the person: "wild type", "JCIM Fig. 4", "redo with Ca²⁺". Tags are
free words, each up to 40 characters, at most 20 a study, the same tag
written once whatever its case; the note is one line of up to 200. Both are
kept in ``study_tags.json`` beside the study's records, never inside them,
so they travel with the study when it is moved, copied or shared, and no
record the software reads is changed. The GUI's All studies page sets them
on a card; an AI app through `fastmdx mcp` reads them and may add tags, and
removing a tag or writing the note stays with the person.

A study the person archived is kept there too (``archived``, when): All
studies leaves it out but in its Archived view. Nothing of the study is
moved or deleted, so its records and the paths in them stay as they were.
"""

from __future__ import annotations

import json
import os
import tempfile
import unicodedata
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

__all__ = ["MOST_TAGS", "NOTE_LENGTH", "TAG_LENGTH", "TAGS_FILE", "add_tags", "set_tags",
           "tags_of"]

TAGS_FILE = "study_tags.json"
TAG_LENGTH = 40
MOST_TAGS = 20
NOTE_LENGTH = 200


def _clean(text: Any, most: int) -> str | None:
    """One line of words: spaces collapsed, no control characters; None
    where it is not text or is too long."""
    if not isinstance(text, str):
        return None
    kept = "".join(" " if unicodedata.category(c).startswith("C") else c for c in text)
    kept = " ".join(kept.split())
    return kept if len(kept) <= most else None


def _kept(folder: str | Path) -> dict[str, Any]:
    try:
        kept = json.loads((Path(folder) / TAGS_FILE).read_text(encoding="utf-8"))
    except (OSError, ValueError):
        kept = {}
    return kept if isinstance(kept, dict) else {}


def _archived_at(kept: dict[str, Any]) -> str | None:
    when = kept.get("archived")
    return when if isinstance(when, str) and when else None


def tags_of(folder: str | Path) -> dict[str, Any]:
    """A study's tags and note, as kept, none where none were given; and
    ``archived`` True where the person archived it."""
    kept = _kept(folder)
    tags = [t for t in (kept.get("tags") or []) if _clean(t, TAG_LENGTH)]
    note = _clean(kept.get("note") or "", NOTE_LENGTH) or ""
    said: dict[str, Any] = {"tags": tags[:MOST_TAGS], "note": note}
    if _archived_at(kept) is not None:
        said["archived"] = True
    return said


def _checked_tags(given: Any) -> tuple[list[str] | None, str | None]:
    if not isinstance(given, list):
        return None, "Tags are given as a list of words."
    tags: list[str] = []
    seen: set[str] = set()
    for tag in given:
        cleaned = _clean(tag, TAG_LENGTH)
        if cleaned is None:
            return None, f"A tag is up to {TAG_LENGTH} characters of text on one line."
        if cleaned and cleaned.casefold() not in seen:
            seen.add(cleaned.casefold())
            tags.append(cleaned)
    if len(tags) > MOST_TAGS:
        return None, f"A study has {MOST_TAGS} tags at most."
    return tags, None


def _write(folder: Path, tags: list[str], note: str,
           archived: bool | None = None) -> dict[str, Any]:
    """Kept, the study archived as it was unless ``archived`` says."""
    from fastmdxplora.gui.browse import is_study

    if not folder.is_dir() or not is_study(folder):
        return {"ok": False, "reason": f"{folder} is not a study."}
    now = datetime.now(timezone.utc).isoformat(timespec="seconds")
    record: dict[str, Any] = {"tags": tags, "note": note, "updated": now}
    when = _archived_at(_kept(folder))
    if archived is not None:
        when = (when or now) if archived else None
    if when:
        record["archived"] = when
    try:
        handle, temporary = tempfile.mkstemp(prefix=".study_tags.", suffix=".tmp", dir=folder)
        with os.fdopen(handle, "w", encoding="utf-8") as out:
            json.dump(record, out, ensure_ascii=False, indent=2)
        os.replace(temporary, folder / TAGS_FILE)
    except OSError as exc:
        return {"ok": False, "reason": f"The study's tags could not be kept in its folder: "
                                       f"{exc.strerror or exc}."}
    said: dict[str, Any] = {"ok": True, "tags": tags, "note": note}
    if when:
        said["archived"] = True
    return said


def set_tags(folder: str | Path, tags: Any, note: Any = None,
             archived: Any = None) -> dict[str, Any]:
    """A study's tags as given, replacing those it had, its note where one
    is given (an empty one clears it), and archived or not where that is
    given (True or False; anything else leaves it as it was)."""
    base = Path(folder)
    checked, reason = _checked_tags(tags)
    if checked is None:
        return {"ok": False, "reason": reason}
    if note is None:
        kept_note = tags_of(base)["note"]
    else:
        kept_note = _clean(note, NOTE_LENGTH)
        if kept_note is None:
            return {"ok": False, "reason": f"A note is one line of up to {NOTE_LENGTH} "
                                           "characters."}
    return _write(base, checked, kept_note, archived if isinstance(archived, bool) else None)


def add_tags(folder: str | Path, tags: Any) -> dict[str, Any]:
    """Tags added to those a study has, the note left as it is."""
    base = Path(folder)
    checked, reason = _checked_tags(tags)
    if checked is None:
        return {"ok": False, "reason": reason}
    now = tags_of(base)
    merged, why = _checked_tags(now["tags"] + checked)
    if merged is None:
        return {"ok": False, "reason": why}
    said = _write(base, merged, now["note"])
    if said.get("ok"):
        known = {t.casefold() for t in now["tags"]}
        said["added"] = [t for t in checked if t.casefold() not in known]
    return said
