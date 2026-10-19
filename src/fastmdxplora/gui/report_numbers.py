"""Each number the report states, traced to the record it came from.

The report is written from the study's records: the analyses' findings,
the Manifest, the resolved Config, the setup and simulation records. On
the Report page a number given with a unit, an error, a percent or what
it counts (atoms, frames, independent samples), and any number with a
decimal point, opens the record that holds it: the file, the place in
it and the value as recorded, converted where the report gives it in
other units (nm and angstrom, fs, ps, ns and microseconds, K and degrees C).

A number is matched to a record value that rounds to it as the report
gives it, in a unit of the same kind. Where several do, the one in the
section's own analysis, or whose name says what the number counts, is
taken. A number no record holds is said so, with the nearest value of
its kind where there is one: a method's own threshold, a citation, a
number the report works out from more than one record, or a report
written before its records last changed. Code alone, no reading of the
prose beyond the words next to the number.

Records are read through `record_access`, so in a read held to the study
a record linked out of it is not read. A number no record that was read
holds is then said as not traced, with the records not read named, never
as one no record holds.
"""

from __future__ import annotations

import fnmatch
import json
import math
import os
import re
from collections.abc import Iterator
from pathlib import Path
from typing import Any

from fastmdxplora import record_access as _records
from fastmdxplora.refusals import OutsideWorkspace

__all__ = ["record_files", "records_of", "traced_html", "trace_number"]

#: The records a report is written from, under the study (or each run of
#: a study of several). Streams and indexes of frames are not records of a
#: result.
_RECORDS = (
    "manifest.json", "stopping.json", "resolved_config.yml",
    "setup/*.json", "simulation/simulation_parameters.json", "simulation/cost.json",
    "analysis/analysis_manifest.json", "analysis/*/options.json",
    "batch_manifest.json", "members.json", "comparison/*.json",
)
#: A list longer than this is a series, not a value the report states.
_LONGEST_LIST = 20
_LARGEST_FILE = 2_000_000

#: Each unit the report writes: its kind, and its size in the kind's base
#: unit (nm, ps, K, bar, kJ/mol, M, ns/day, 1/ps, nm^2, degree, amu nm^2).
_UNITS: dict[str, tuple[str, float]] = {
    "amu nm\u00b2": ("inertia", 1.0), "amu nm^2": ("inertia", 1.0),
    "nm\u00b2": ("area", 1.0), "nm^2": ("area", 1.0), "\u00c5\u00b2": ("area", 0.01),
    "ns/day": ("rate", 1.0),
    "kJ/mol": ("energy", 1.0), "kJ mol\u207b\u00b9": ("energy", 1.0),
    "kcal/mol": ("energy", 4.184), "kcal mol\u207b\u00b9": ("energy", 4.184),
    "ps\u207b\u00b9": ("per_time", 1.0), "ps^-1": ("per_time", 1.0), "/ps": ("per_time", 1.0),
    "nm": ("length", 1.0), "\u00c5": ("length", 0.1),
    "fs": ("time", 0.001), "ps": ("time", 1.0), "ns": ("time", 1000.0),
    "\u00b5s": ("time", 1.0e6), "\u03bcs": ("time", 1.0e6),
    "\u00b0C": ("temperature", 1.0), "K": ("temperature", 1.0),
    "bar": ("pressure", 1.0), "atm": ("pressure", 1.01325),
    "mM": ("concentration", 0.001), "M": ("concentration", 1.0),
    "%": ("percent", 1.0), "\u00b0": ("angle", 1.0),
}
#: A record's unit as the end of its key names it.
_SUFFIXES: tuple[tuple[str, str], ...] = (
    ("_per_ps", "ps\u207b\u00b9"), ("ns_per_day", "ns/day"), ("_kjmol", "kJ/mol"),
    ("_kcalmol", "kcal/mol"), ("_nm2", "nm\u00b2"), ("_nm", "nm"), ("_angstrom", "\u00c5"),
    ("_A", "\u00c5"), ("_fs", "fs"), ("_ps", "ps"), ("_ns", "ns"), ("_us", "\u00b5s"),
    ("_K", "K"), ("_kelvin", "K"), ("_bar", "bar"), ("_atm", "atm"), ("_M", "M"),
    ("_mM", "mM"), ("_deg", "\u00b0"), ("_degrees", "\u00b0"),
)
#: A record's unit as a word in its key says it.
_WORDS: tuple[tuple[str, str], ...] = (
    ("temperature", "K"), ("pressure", "bar"), ("ns_per_day", "ns/day"),
)
#: What a unit a record gives in a field of its own is called here.
_NAMED = {"A": "\u00c5", "angstrom": "\u00c5", "nm^2": "nm\u00b2", "A^2": "\u00c5\u00b2",
          "kj/mol": "kJ/mol", "kcal/mol": "kcal/mol", "deg": "\u00b0", "degrees": "\u00b0",
          "us": "\u00b5s", "1/ps": "ps\u207b\u00b9", "amu*nm^2": "amu nm\u00b2"}
#: The fields of a record's dict that the dict's own `unit` applies to.
_VALUED = {"mean", "value", "standard_error", "standard_deviation", "error", "sem", "std",
           "median", "minimum", "maximum", "min", "max", "lower", "upper", "first", "last"}
#: What a count's noun is called in a record's keys.
_COUNTED = {
    "atoms": ("atom",), "frames": ("frame", "discard"), "residues": ("residue",),
    "samples": ("sample", "independent"), "replicas": ("replica", "member", "run"),
    "steps": ("step",), "bonds": ("bond",), "clusters": ("cluster",), "ions": ("ion",),
    "chains": ("chain",), "molecules": ("molecule",), "waters": ("water",),
    "runs": ("run", "member"), "contacts": ("contact",), "pieces": ("segment", "piece"),
}

_NUMBER = re.compile(
    r"(?<![\w.,/:#\-\u2212])([-\u2212]?)(\d{1,3}(?:,\d{3})+|\d+)(\.\d+)?(?:[eE]([-+]?\d+))?"
    r"(?![\w\u00b2\u00b3]|[.,]\d|/\w)")
_UNIT = re.compile(r"\s?(" + "|".join(re.escape(u) for u in sorted(_UNITS, key=len, reverse=True))
                   + r")(?![\w\u00b2\u00b3])")
_ERROR = re.compile(r"\s*(?:\u00b1|\+/-|&plusmn;)\s*")
_NOUN = re.compile(r"\s+(?:[a-z][a-z-]*\s+){0,2}?(" + "|".join(_COUNTED) + r")\b")
#: Inside these the text is not the report's statement of a number.
_SKIPPED = {"code", "pre", "a", "h1", "h2", "h3", "h4", "h5", "h6", "script", "style", "sup"}


# ---------------------------------------------------------------------------
# The records
# ---------------------------------------------------------------------------

def _unit_of_key(key: str) -> str:
    for suffix, unit in _SUFFIXES:
        if key.endswith(suffix):
            return unit
    lowered = key.lower()
    for word, unit in _WORDS:
        if word in lowered:
            return unit
    return ""


def _named_unit(said: Any) -> str:
    if not isinstance(said, str) or not said.strip():
        return ""
    said = said.strip()
    if said in _UNITS:
        return said
    return _NAMED.get(said, _NAMED.get(said.lower(), ""))


def _leaves(value: Any, path: str, unit: str) -> Iterator[tuple[str, float, str]]:
    if isinstance(value, bool):
        return
    if isinstance(value, (int, float)):
        if math.isfinite(value):
            yield path, float(value), unit
        return
    if isinstance(value, dict):
        own = _named_unit(value.get("unit") or value.get("units"))
        for key, item in value.items():
            key = str(key)
            here = _unit_of_key(key) or (own if key in _VALUED else "")
            yield from _leaves(item, f"{path}.{key}" if path else key, here)
        return
    if isinstance(value, list) and len(value) <= _LONGEST_LIST:
        for i, item in enumerate(value):
            yield from _leaves(item, f"{path}[{i}]", unit)


def _read(path: Path, base: Path, not_read: list[str] | None) -> Any:
    try:
        if _records.stat(path).st_size > _LARGEST_FILE:
            return None
        text = _records.read_text(path, encoding="utf-8")
        if path.suffix in (".yml", ".yaml"):
            import yaml

            return yaml.safe_load(text)
        return json.loads(text)
    except OutsideWorkspace:
        # Refused in a read held to the study: not read, which is not a
        # record that holds nothing.
        if not_read is not None:
            not_read.append(path.relative_to(base).as_posix())
        return None
    except (OSError, ValueError, ImportError):
        return None
    except Exception:  # noqa: BLE001 - a YAML error is not an OSError; the record is skipped
        return None


def records_of(root: Path | str, not_read: list[str] | None = None) -> list[dict[str, Any]]:
    """Every number the study's records hold: its file (from the study's
    folder), its place in the file, its value and its unit where the record
    says it. A study of several runs gives its own and each run's. Each
    record refused (`record_access`) is named in ``not_read``."""
    base = Path(root)
    found: list[dict[str, Any]] = []
    for path in record_files(base, not_read):
        data = _read(path, base, not_read)
        if data is None:
            continue
        file = path.relative_to(base).as_posix()
        for key, value, unit in _leaves(data, "", ""):
            found.append({"file": file, "key": key, "value": value, "unit": unit})
    return found


# ---------------------------------------------------------------------------
# One number
# ---------------------------------------------------------------------------

def _in_base(value: float, unit: str) -> tuple[str, float] | None:
    if unit not in _UNITS:
        return None
    kind, size = _UNITS[unit]
    if unit == "\u00b0C":
        return kind, value + 273.15
    return kind, value * size


def _from_base(kind: str, value: float, unit: str) -> float:
    if unit == "\u00b0C":
        return value - 273.15
    return value / _UNITS[unit][1]


def _figures(text: str) -> int:
    mantissa = re.split(r"[eE]", text)[0]
    digits = mantissa.replace(",", "").replace(".", "").lstrip("-\u2212").lstrip("0")
    return len(digits)


def _decimals(text: str) -> int:
    """The place a number as written is given to: 0.0122 to the fourth
    decimal, 9.5e-05 to the sixth, 7,439 to the unit."""
    mantissa, _, exponent = text.lower().partition("e")
    places = len(mantissa.split(".")[1]) if "." in mantissa else 0
    return places - (int(exponent) if exponent else 0)


def _rounds_to(value: float, said: str, decimals: int) -> tuple[bool, bool]:
    """Whether ``value`` rounds to ``said`` as given, and whether it is
    the number the software itself prints for it (four figures, its
    trailing zeros dropped): 0.02 is 0.020003 so printed, and 0.02089
    only rounded."""
    from fastmdxplora.statistics import four_figures

    number = float(said.replace(",", "").replace("\u2212", "-"))
    step = 0.5 * 10 ** (-decimals)
    near = abs(value - number) <= step * (1 + 1e-9) + 1e-12
    if not near:
        return False, False
    plain = said.replace("\u2212", "-")
    return True, four_figures(value) == plain


#: Words next to a number that say nothing of which record it is.
_PLAIN_WORDS = {
    "the", "and", "for", "with", "from", "was", "were", "are", "its", "this", "that", "each",
    "over", "after", "before", "under", "into", "than", "then", "run", "one", "all", "any",
    "not", "but", "per", "has", "had", "have", "about", "which", "where", "when", "part",
}
#: A record of what is missing or asked for, taken for a number only where
#: the words before it say so.
_ASKED = ("shortfall", "target", "more_")
_ASKING = ("more", "at least", "target")


#: A word for a record's field as a table's head or the prose says it.
_SAYS = {"independent": "effective", "discarded": "discard", "uncertainty": "error"}


def _stem(word: str) -> str:
    word = _SAYS.get(word, word)
    for ending in ("ing", "ed", "es", "s"):
        if word.endswith(ending) and len(word) - len(ending) >= 3:
            return word[:-len(ending)]
    return word


def _words(text: str) -> set[str]:
    found = {w for w in re.findall(r"[a-z][a-z0-9]+", text.lower()) if w not in _PLAIN_WORDS}
    return {_stem(w) for w in found if len(w) >= 3 or w == "ph"}


def _key_words(key: str) -> set[str]:
    return {_stem(w) for w in re.split(r"[._\[\]]+", key.lower()) if w}


def _about(subject: str, record: dict[str, Any]) -> bool:
    file = "/" + record["file"].lower()
    return f"/{subject}/" in file or f".{subject}." in "." + record["key"].lower() + "."


class _Study:
    """A study's records, and what its headings and table rows can name:
    each analysis by its folder (or its title), each run by its folder,
    and each quantity a record is kept under (`potential_energy`)."""

    def __init__(self, records: list[dict[str, Any]]) -> None:
        self.records = records
        self.analyses: set[str] = set()
        self.runs: set[str] = set()
        self.quantities: set[str] = set()
        self.vocabulary: set[str] = set()
        for record in records:
            parts = record["file"].split("/")
            if parts[0] == "runs" and len(parts) > 2:
                self.runs.add(parts[1].lower())
            if "analysis" in parts and parts.index("analysis") + 2 < len(parts):
                self.analyses.add(parts[parts.index("analysis") + 1].lower())
            self.quantities.update(w for w in record["key"].lower().split(".")
                                   if not w.endswith("]"))
            self.vocabulary |= _key_words(record["key"])
        self.titles = _titles(self.analyses)

    def named(self, words: str) -> str:
        said = _named(words)
        return said if said in self.analyses | self.runs | self.quantities else ""

    def titled(self, words: str) -> list[str]:
        """The analyses a column's head names, by folder or by title."""
        plain = " " + re.sub(r"[^a-z0-9_]+", " ", words.lower()) + " "
        return [name for name, titles in self.titles.items()
                if f" {name} " in plain or any(f" {title} " in plain for title in titles)]


def _titles(analyses: set[str]) -> dict[str, list[str]]:
    """Each analysis by the titles the pages give it."""
    titles: dict[str, list[str]] = {name: [] for name in analyses}
    try:
        from fastmdxplora.gui.report_dashboard import ANALYSIS_THEMES
        from fastmdxplora.gui.series import SERIES

        said = [pair for _theme, pairs in ANALYSIS_THEMES for pair in pairs]
        said += [(name, title) for name, (title, _over) in SERIES.items()]
    except Exception:  # noqa: BLE001 - the folders' names alone still serve
        said = []
    for name, title in said:
        if name in titles:
            titles[name].append(re.sub(r"[^a-z0-9_]+", " ", title.lower()).strip())
    return titles


def records_said(study: Any) -> _Study:
    return study if isinstance(study, _Study) else _Study(study)


def trace_number(said: str, unit: str, records: Any, *, counts: str = "",
                 subjects: tuple[str, ...] = (), context: str = "", head: str = "",
                 error: bool = False) -> dict[str, Any]:
    """The record a number as the report gives it comes from, or the
    nearest of its kind where none does.

    ``said`` is the number as written ("7,439", "0.07415"), ``unit`` the
    unit after it, ``counts`` the noun it counts, ``subjects`` the
    analyses, runs or quantities its section, table row or column is
    about (another analysis's or run's records are then not its),
    ``context`` the words before it, ``head`` its table column's head (a
    record must then be what the head names, and the head's unit is the
    number's where the cell gives none), ``error`` whether it is an error
    after a plus-minus."""
    study = records_said(records)
    decimals = _decimals(said)
    figures = _figures(said)
    number = float(said.replace(",", "").replace("−", "-"))
    stems = _COUNTED.get(counts, ())
    if not unit and head:
        in_head = re.search(r"\(([^)]+)\)\s*$", head)
        unit = in_head.group(1).strip() if in_head else ""
    said_unit = unit if unit in _UNITS else ""
    near_words = _words(context + " " + head)
    head_words = _words(head) & study.vocabulary
    asking = any(word in context.lower() for word in _ASKING)
    best: tuple[float, dict[str, Any]] | None = None
    others = 0
    nearest: tuple[tuple[int, float], dict[str, Any]] | None = None
    for record in study.records:
        key = record["key"].lower()
        file = record["file"].lower()
        in_analysis = "analysis/" in file and file.count("/") >= 2
        about = 0
        apart = False
        for subject in subjects:
            if subject in study.runs:
                if file.startswith("runs/"):
                    if file.startswith(f"runs/{subject}/"):
                        about += 1
                    else:
                        apart = True
            elif _about(subject, record):
                about += 1
            elif in_analysis:
                apart = True
        if apart:
            continue
        if not asking and any(word in key for word in _ASKED):
            continue
        key_words = _key_words(key)
        if head_words and not head_words & key_words:
            continue
        relevant = 3 * about
        if stems and any(stem in key for stem in stems):
            relevant += 3
        if error and "error" in key:
            relevant += 2
        named = bool(near_words & key_words)
        if named:
            relevant += 2
        record_unit = record["unit"]
        if said_unit == "%":
            if record_unit:
                continue
            candidates = [record["value"], record["value"] * 100.0]
        elif said_unit:
            if record_unit:
                theirs = _in_base(record["value"], record_unit)
                if theirs is None or theirs[0] != _UNITS[said_unit][0]:
                    continue
                candidates = [_from_base(theirs[0], theirs[1], said_unit)]
            else:
                candidates = [record["value"]]
        else:
            # A number said without its unit (a table's column of means)
            # is a record's with one only where it is plainly that record.
            if record_unit and (counts or not (about and named)):
                continue
            candidates = [record["value"]]
        for value in candidates:
            near, exact = _rounds_to(value, said, decimals)
            if not near:
                if relevant >= 2 and number and (record_unit or not said_unit):
                    gap = abs(value - number) / abs(number)
                    rank = (-relevant, gap)
                    if gap < 0.5 and (nearest is None or rank < nearest[0]):
                        nearest = (rank, {**record, "as_said": value})
                continue
            # A unit the record does not say, or a number with no unit or
            # noun: taken only where the record's name says what it is, or
            # the number is too precise to meet by chance.
            known = bool(record_unit) and bool(said_unit) and said_unit != "%"
            if not known and not relevant and not (exact and figures >= 3):
                continue
            # Of two alike, the shorter name: `findings.mean.n_frames`
            # before `findings.hydrophobic_sasa.n_frames`.
            score = 10 * exact + 5 * known + relevant - len(key) / 1000.0
            if best is None or score > best[0]:
                if best is not None:
                    others += 1
                best = (score, {**record, "as_said": value})
            else:
                others += 1
    if best is not None:
        return {"traced": True, "record": _said_record(best[1], said_unit), "others": others}
    return {"traced": False, "nearest": _said_record(nearest[1], said_unit) if nearest else None}


def _said_record(record: dict[str, Any], unit: str) -> dict[str, Any]:
    value = record["value"]
    said = repr(value) if not float(value).is_integer() else str(int(value))
    return {"file": record["file"], "key": record["key"], "value": said,
            "unit": record["unit"], "in": unit, "as_said": record["as_said"]}


# ---------------------------------------------------------------------------
# The report's HTML
# ---------------------------------------------------------------------------

def _named(words: str) -> str:
    """What a heading or a row names, as a folder or a key is named:
    "End_To_End" is end_to_end, "potential energy" potential_energy."""
    plain = re.sub(r"<[^>]+>", "", words).replace("\\_", "_").strip().lower()
    return re.sub(r"[^a-z0-9_]+", "_", plain).strip("_")


def _numbers_in(text: str, every: bool) -> Iterator[tuple[int, int, str, str, str, bool]]:
    """Each number in a run of text the report states as a quantity: its
    span, as written, its unit, what it counts, and whether it is an
    error. ``every``: a table's cell, where any number is one."""
    previous_end = -1
    for found in _NUMBER.finditer(text):
        start, end = found.span()
        said = text[start:end]
        before = text[max(0, previous_end):start]
        error = bool(_ERROR.fullmatch(before)) if previous_end >= 0 else False
        previous_end = end
        rest = text[end:]
        after_error = _ERROR.match(rest)
        has_error = False
        if after_error:
            next_number = _NUMBER.match(rest, after_error.end())
            if next_number:
                has_error = True
                rest = rest[next_number.end():]
        unit_found = _UNIT.match(rest)
        unit = unit_found.group(1) if unit_found else ""
        noun_found = None if unit else _NOUN.match(text[end:])
        counts = noun_found.group(1) if noun_found else ""
        if not (every or unit or counts or "." in said or error or has_error):
            continue
        yield start, end, said, unit, counts, error


def traced_html(html: str, root: Path | str,
                not_read: list[str] | None = None) -> tuple[str, list[dict[str, Any]]]:
    """The report's HTML with each number it states marked, and what each
    traced to. A mark is `<span class="report-number" data-number="i">`
    around the number alone. A number not traced where a record was not
    read names those records (``not_read``, filled here too): it may be
    one of theirs."""
    refused: list[str] = []
    study = _Study(records_of(root, refused))
    unread = sorted(set(refused))
    if not_read is not None:
        not_read.extend(unread)
    parts = re.split(r"(<[^>]+>)", html)
    out: list[str] = []
    numbers: list[dict[str, Any]] = []
    open_tags: list[str] = []
    section = ""
    heading: list[str] | None = None
    # A table's column heads, and the row's first cell: what its numbers are.
    heads: list[str] = []
    row = ""
    cell: list[str] | None = None
    column = -1
    in_head = False
    # A paragraph or list item led by a name in bold ("**rg.** The
    # replicas spread ...") is about what it names.
    lead = ""
    bold: list[str] | None = None
    for part in parts:
        if part.startswith("<"):
            name = re.match(r"</?\s*([a-zA-Z0-9]+)", part)
            tag = name.group(1).lower() if name else ""
            closing = part.startswith("</")
            if tag in ("p", "li"):
                lead = ""
            elif tag == "strong":
                if closing and bold is not None:
                    named = study.named("".join(bold))
                    lead = lead or named
                    bold = None
                elif not closing:
                    bold = []
            if tag == "table" and not closing:
                heads = []
            elif tag == "thead":
                in_head = not closing
            elif tag == "tr" and not closing:
                row, column = "", -1
            elif tag in ("td", "th"):
                if closing:
                    if cell is not None:
                        words = re.sub(r"<[^>]+>", "", "".join(cell)).strip()
                        if in_head or tag == "th":
                            heads.append(words)
                        elif column == 0:
                            row = study.named(words)
                    cell = None
                else:
                    column += 1
                    cell = []
            if closing:
                if tag in open_tags:
                    while open_tags and open_tags.pop() != tag:
                        pass
                if tag in ("h2", "h3") and heading is not None:
                    # A section of the run (Methods, Results) is no analysis.
                    said = study.named("".join(heading))
                    section = said if tag == "h3" or said in study.analyses else ""
                    heading = None
            elif not part.endswith("/>") and tag in _SKIPPED:
                open_tags.append(tag)
                if tag in ("h2", "h3"):
                    heading = []
            out.append(part)
            continue
        if heading is not None:
            heading.append(part)
        if cell is not None:
            cell.append(part)
        if bold is not None:
            bold.append(part)
        if open_tags or not part or (cell is not None and column == 0):
            out.append(part)
            continue
        in_cell = cell is not None
        head = heads[column] if in_cell and 0 <= column < len(heads) else ""
        subjects = tuple(s for s in ([row] if in_cell and row else [lead or section])
                         + (study.titled(head) if head else []) if s)
        cursor = 0
        pieces: list[str] = []
        for start, end, said, unit, counts, error in _numbers_in(part, in_cell):
            trace = trace_number(said, unit, study, counts=counts, subjects=subjects,
                                 context=part[max(0, start - 48):start], head=head,
                                 error=error)
            index = len(numbers)
            if unread and not trace["traced"]:
                trace = {**trace, "not_read": unread}
            numbers.append({"said": said + (" " + unit if unit else ""), "counts": counts,
                            "subjects": list(subjects), **trace})
            pieces.append(part[cursor:start])
            pieces.append(f'<span class="report-number" data-number="{index}" '
                          f'data-traced="{str(trace["traced"]).lower()}" role="button" '
                          f'tabindex="0" aria-expanded="false">{said}</span>')
            cursor = end
        pieces.append(part[cursor:])
        out.append("".join(pieces))
    return "".join(out), numbers


def record_files(root: Path | str, not_read: list[str] | None = None) -> list[Path]:
    """The record files a report is traced against, for knowing when they
    changed. Listed through `record_access`: in a read held to the study,
    one there that is linked out of it (or a folder of them) is not listed
    but named in ``not_read``."""
    base = Path(root)
    runs = [p for p in _matching(base, "runs/*", base, not_read, folders=True)]
    found: list[Path] = []
    for folder in [base] + runs:
        for pattern in _RECORDS:
            found.extend(_matching(folder, pattern, base, not_read))
    return found


def _matching(folder: Path, pattern: str, base: Path, not_read: list[str] | None, *,
              folders: bool = False) -> list[Path]:
    """The files (or ``folders``) under ``folder`` a pattern of names
    matches, each part a name or a wildcard; one that is there but may not
    be read named in ``not_read``."""
    def refused(path: Path) -> bool:
        if os.path.lexists(path) and not _records.can_read(path):
            if not_read is not None:
                not_read.append(path.relative_to(base).as_posix())
            return True
        return False

    places = [folder]
    parts = pattern.split("/")
    for index, part in enumerate(parts):
        last = index == len(parts) - 1
        found: list[Path] = []
        for place in places:
            if any(char in part for char in "*?["):
                try:
                    children = sorted(child for child in _records.iterdir(place)
                                      if fnmatch.fnmatchcase(child.name, part))
                except OSError:
                    children = []
            else:
                children = [place / part]
            for child in children:
                wanted = _records.is_file if last and not folders else _records.is_dir
                if wanted(child):
                    found.append(child)
                else:
                    refused(child)
        places = found
    return places
