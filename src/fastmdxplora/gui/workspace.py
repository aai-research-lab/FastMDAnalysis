"""The studies in a workspace, as cards, and two of them compared.

The GUI showed one study at a time, and another was opened by walking to
its folder. This finds the studies under a folder and says each in a card:
its structure, what kind of study it is, where it stands, when it was made,
and the means it recorded with their errors, with a figure it plotted. Two
chosen are compared as `fastmdx diff` compares them, the settings in which
they differ with defaults filled in, and beside that the means both
recorded, a difference marked only past twice their combined error, the
rule the Analysis page and the written comparison use.

Everything is read from the studies' own records; nothing is written.
"""

from __future__ import annotations

import json
import math
import os
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from fastmdxplora.system_id import system_id

__all__ = ["studies_in", "card_of", "head_of", "studies_compared", "compared_side", "sides_compared",
           "means_of", "tags_used", "thumbnail_of", "thumbnail_series", "study_key", "study_keyed",
           "study_state"]

#: A difference is marked past this many combined standard errors.
RESOLVED_AT = 2.0

#: How deep under the folder to look, and how many studies to say: a home
#: folder can hold thousands of folders, and the page has to answer.
DEEPEST = 3
MOST = 200

#: Folders a study keeps inside it that are never studies of their own.
_INSIDE_A_STUDY = {"runs", "superseded", "setup", "simulation", "analysis", "report",
                   "comparison", "shared_setup", "free_energy", "pmf"}

#: How many folders a study named by its key is looked for in.
_LOOKED_IN = 20000

#: The measures a card names first, in order.
_FIRST = ("rmsd", "rg", "sasa", "rmsf", "hbonds", "ligand_rmsd")


def studies_in(root: Path | str, *, deepest: int = DEEPEST,
               most: int = MOST) -> dict[str, Any]:
    """The studies under ``root``, newest first, each as a card."""
    from fastmdxplora.gui.browse import is_study

    base = Path(root).expanduser()
    try:
        base = base.resolve()
    except OSError:
        pass
    if not base.is_dir():
        return {"ok": False, "reason": f"No such folder: {base}"}
    found: list[Path] = []
    more = False

    def walk(folder: Path, depth: int) -> None:
        nonlocal more
        if len(found) >= most:
            more = True
            return
        try:
            children = sorted(p for p in folder.iterdir()
                              if p.is_dir() and not p.name.startswith("."))
        except OSError:
            return
        for child in children:
            if len(found) >= most:
                more = True
                return
            if child.name in _INSIDE_A_STUDY:
                continue
            if is_study(child):
                found.append(child)
            elif depth < deepest:
                walk(child, depth + 1)

    if is_study(base):
        found.append(base)
    else:
        walk(base, 1)
    cards = [card_of(folder) for folder in found]
    cards.sort(key=lambda card: card.get("when") or "", reverse=True)
    return {"ok": True, "root": str(base), "studies": cards, "more": more,
            "tags_used": tags_used(cards)}


def study_key(folder: Path | str) -> str:
    """The key a page names a study by in its address: the study's folder,
    hashed, so no path is said; the same key `/api/status` gives."""
    import hashlib

    return hashlib.sha256(str(Path(folder).expanduser().resolve()).encode("utf-8")).hexdigest()[:16]


def study_keyed(key: str, near: list[Path | str | None],
                under: list[Path | str | None] = ()) -> Path | None:
    """The study an address names by its key: one of the folders ``near``
    (the one open, the one running), the study each sits in, or a run
    inside one of those; else a study under the folders ``under`` (the ones
    the GUI looks in) as All studies finds them, or a run inside one. None
    where none has that key. Only ``under`` is walked: the folders around a
    study are not looked through."""
    import re

    from fastmdxplora.gui.browse import is_study

    if not isinstance(key, str) or not re.fullmatch(r"[0-9a-f]{16}", key):
        return None
    looked = 0
    seen: set[Path] = set()

    def matches(folder: Path) -> bool:
        if folder in seen:
            return False
        seen.add(folder)
        return is_study(folder) and study_key(folder) == key

    def runs_of(study: Path) -> Path | None:
        try:
            runs = [p for p in (study / "runs").iterdir() if p.is_dir()]
        except OSError:
            return None
        return next((run for run in runs if matches(run)), None)

    def walk(folder: Path, depth: int) -> Path | None:
        nonlocal looked
        try:
            children = sorted(p for p in folder.iterdir()
                              if p.is_dir() and not p.name.startswith("."))
        except OSError:
            return None
        for child in children:
            looked += 1
            if looked > _LOOKED_IN:
                return None
            if child.name in _INSIDE_A_STUDY:
                continue
            if matches(child):
                return child
            if is_study(child):
                found = runs_of(child)
            elif depth < DEEPEST:
                found = walk(child, depth + 1)
            else:
                found = None
            if found is not None:
                return found
        return None

    def resolved(folders: Any) -> list[Path]:
        kept: list[Path] = []
        for given in folders:
            if given is None:
                continue
            try:
                kept.append(Path(given).expanduser().resolve())
            except OSError:
                continue
        return kept

    bases: list[Path] = []
    for base in resolved(near) + resolved(under):
        # A run inside a study: the study, and its other runs, are near too.
        bases.extend(b for b in (base, base.parent.parent) if b not in bases)
    for base in bases:
        if not base.is_dir():
            continue
        if matches(base):
            return base
        found = runs_of(base) if is_study(base) else None
        if found is not None:
            return found
    for base in resolved(under):
        if base.is_dir():
            found = walk(base, 1)
            if found is not None:
                return found
    return None


def tags_used(cards: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """The tags the studies carry, the most used first, each with how many
    carry it: one tag whatever its case, written as it was first."""
    counted: dict[str, list[Any]] = {}
    for card in cards:
        for tag in card.get("tags") or []:
            counted.setdefault(tag.casefold(), [tag, 0])[1] += 1
    return [{"tag": tag, "studies": n}
            for tag, n in sorted(counted.values(), key=lambda kept: (-kept[1], kept[0].casefold()))]


def head_of(folder: Path | str) -> dict[str, Any]:
    """What a study is, from its records: its folder, its system, what kind
    of study it is and where it stands. Read through `record_access`, held
    to the study where the read is (the Compare page): where the record its
    state is said from was not read, the state is not known, never one
    taken from that record's absence."""
    from fastmdxplora import record_access

    base = Path(folder)
    head = _head(base, *_records_of(base))
    for name in ("batch_manifest.json", "manifest.json"):
        if os.path.lexists(base / name) and not record_access.can_read(base / name):
            head["state"] = f"not known: {name} was not read"
            break
    return head


def _head(base: Path, config: Any, batch: Any, manifest: Any) -> dict[str, Any]:
    named = _system_of(config, batch, manifest)
    return {
        "path": str(base),
        "name": base.name,
        # The system in four capital characters (fastmdxplora.system_id),
        # and what named it: a PDB entry, or its structure file's name.
        "system": system_id(named),
        "structure": named,
        "kind": _kind_of(config, batch, simulated=_simulated(base)),
        "state": _state_of(base, batch, manifest),
    }


def _records_of(base: Path) -> tuple[Any, Any, Any]:
    """The study's Config (its phase settings lifted), batch manifest and
    Manifest, each None where it has none."""
    config = _lifted(_read_yaml(base / "resolved_config.yml")
                     or _read_yaml(base / "exploration.yml"))
    return config, _read_json(base / "batch_manifest.json"), _read_json(base / "manifest.json")


def card_of(folder: Path | str) -> dict[str, Any]:
    """What a card says of one study, from its records."""
    base = Path(folder)
    config, batch, manifest = _records_of(base)
    card: dict[str, Any] = {
        **_head(base, config, batch, manifest),
        "when": _when(base, manifest),
        # When anything was last done with it: run, analysed, reported,
        # tagged, viewed or talked about.
        "last_active": _last_active(base),
        "means": [],
        # A figure it plotted, else a picture of its backbone, else none.
        "thumbnail": _picture_of(base),
        # The measure whose series the card plots in the page's colours,
        # where it has one: its figure is white on any scheme.
        "series": _series_of(base),
    }
    from fastmdxplora.study_tags import tags_of

    # What the person said the study is, kept beside its records.
    card.update(tags_of(base))
    simulation = (config or {}).get("simulation") if isinstance(config, dict) else None
    from fastmdxplora.simulation.resume import extended_production

    # An extended study's config carries its last piece's length after the
    # join, so its pieces are what say how long it ran.
    extended = extended_production(base) if not isinstance(batch, dict) else None
    if extended:
        card["production_ns"], card["pieces"] = round(extended[0], 6), extended[1]
    elif isinstance(simulation, dict) and simulation.get("duration_ns") is not None:
        card["production_ns"] = simulation.get("duration_ns")
    forcefield = ((config or {}).get("setup") or {}).get("forcefield") \
        if isinstance(config, dict) else None
    if forcefield and str(forcefield) != "auto":
        card["forcefield"] = str(forcefield)
    if not isinstance(batch, dict):
        card["means"] = [{"analysis": name, **mean} for name, mean in _means(base)[:3]]
    pmf = _read_json(base / "pmf.json")
    if isinstance(pmf, dict):
        # Whether the windows gave a free energy, and if not the first
        # sentence of why: the profile itself is on the study's pages.
        refused = str(pmf.get("refused") or "")
        card["free_energy"] = ({"refused": refused.split(". ")[0].rstrip(".") + "."}
                               if refused else {"recombined": True})
    return card


def _picture_of(base: Path) -> str | None:
    if thumbnail_of(base) is not None:
        return "figure"
    from fastmdxplora.gui.backbone_picture import structure_for_picture

    return "backbone" if structure_for_picture(base) is not None else None


def _series_of(base: Path) -> str | None:
    """The first of the card's measures with a series to plot (series.py)."""
    from fastmdxplora.gui.series import SERIES

    for name in _FIRST:
        if name in SERIES and (base / "analysis" / name / f"{name}.dat").is_file():
            return name
    return None


#: The most points a card's series is sent: a card is 300 px across.
CARD_POINTS = 240


def thumbnail_series(folder: Path | str) -> dict[str, Any]:
    """The series a card plots, thinned to what a card can show."""
    from fastmdxplora.gui.series import series_payload

    base = Path(folder)
    name = _series_of(base)
    if name is None:
        return {"ok": False, "reason": "no series to plot"}
    found = series_payload(base, name)
    if not found.get("ok"):
        return found
    n = len(found["y"])
    every = max(1, -(-n // CARD_POINTS))
    mean = found.get("mean") or {}
    return {
        "ok": True, "analysis": name, "label": found["label"], "unit": found["unit"],
        "kind": found["kind"], "x_label": found.get("x_label") or "",
        "x": found["x"][::every], "y": found["y"][::every],
        "mean": mean.get("value"), "from_x": mean.get("from_x"),
    }


def thumbnail_of(folder: Path | str) -> Path | None:
    """A figure the study plotted, for its card: the first of its measures'
    own figures, else the report's summary figure."""
    base = Path(folder)
    for name in _FIRST:
        figure = base / "analysis" / name / f"{name}.png"
        if figure.is_file():
            return figure
    summary = base / "report" / "analysis_summary.png"
    if summary.is_file():
        return summary
    # An umbrella study's own figure is its free energy profile.
    for profile in sorted((base / "pmf").glob("*.png")) if (base / "pmf").is_dir() else []:
        return profile
    analysis = base / "analysis"
    if analysis.is_dir():
        for figure in sorted(analysis.glob("*/*.png")):
            return figure
    return None


def studies_compared(first: Path | str, second: Path | str) -> dict[str, Any]:
    """Two studies side by side: the settings in which they differ, and the
    means both recorded, a difference marked only where it is resolved."""
    a, b = Path(first), Path(second)
    try:
        one, other = compared_side(a), compared_side(b)
    except Exception as exc:  # noqa: BLE001 - said, not raised
        from fastmdxplora.refusals import refusal_of

        return {"ok": False, "reason": refusal_of(exc).message}
    return {"ok": True, "first": card_of(a), "second": card_of(b),
            **sides_compared(one, other)}


def compared_side(base: Path) -> dict[str, Any]:
    """What a comparison reads of one study: its Config, the phases and
    analyses it ran, and the means it recorded. Each record is read through
    `record_access`, so a read held to the study refuses one linked out of
    it; an analysis whose record was refused is named in ``not_read``,
    never taken for one that recorded no mean. Raises where the study has
    no Config to compare."""
    from fastmdxplora import record_access
    from fastmdxplora.config.diff import config_of

    config, _ = config_of(base)
    not_read: list[str] = []
    try:
        folders = [f for f in record_access.iterdir(base / "analysis") if record_access.is_dir(f)]
    except OSError:
        folders = []
    for folder in folders:
        options = folder / "options.json"
        if os.path.lexists(options) and not record_access.can_read(options):
            not_read.append(folder.name)
    return {"config": config, "ran": _phases_ran(base), "analysed": _analysed(base),
            "ran_analyses": _ran_analyses(base), "means": means_of(base),
            "not_read": sorted(not_read)}


def means_of(base: Path) -> dict[str, dict[str, Any]]:
    """Every mean a comparison sets beside another's: each analysis's own,
    then the thermodynamic means as the run's Overview reads them (from the
    thermodynamics analysis where it recorded them, from production's
    energy file otherwise), each with its error only where it is
    determined. A thermodynamic mean whose record was not read is said so."""
    means = dict(_means(base))
    try:
        from fastmdxplora.gui.overview_view import THERMODYNAMICS, _thermodynamics

        found = _thermodynamics(base)
    except Exception:  # noqa: BLE001 - compared without them
        return means
    shown = found.get("means") if isinstance(found.get("means"), dict) else {}
    for key, label, unit in THERMODYNAMICS:
        said = shown.get(key)
        if isinstance(said, dict) and _finite(said.get("mean")) is not None:
            determined = bool(said.get("determined"))
            means[key] = {"label": label, "mean": _finite(said.get("mean")),
                          "error": _finite(said.get("error")) if determined else None,
                          "unit": unit,
                          "withheld": None if determined else (str(said.get("why") or "")
                                                               or "not determined")}
        elif key in shown and found.get("not_read"):
            means[key] = {"label": label, "mean": None, "error": None, "unit": unit,
                          "withheld": None, "not_read": True}
        elif key in (found.get("no_mean") or {}):
            # Its analysis ran and said why it gave no mean: not determined,
            # never "not recorded".
            means[key] = {"label": label, "mean": None, "error": None, "unit": unit,
                          "withheld": found["no_mean"][key]}
    return means


def sides_compared(one: dict[str, Any], other: dict[str, Any],
                   errors: dict[str, tuple[float, float]] | None = None) -> dict[str, Any]:
    """The settings in which two studies differ and the means each recorded,
    from what `compared_side` read of each. ``errors`` gives, by quantity,
    the error of each one's mean to compare them with in place of the
    error each run estimated for itself: two replicas' (the Compare page,
    `batch.aggregate.replica_error`)."""
    from fastmdxplora.config.diff import differences

    found = [d for d in differences(one["config"], other["config"],
                                    ran=(one["ran"], other["ran"]))
             if not d.where_only]

    def not_reached(d: Any) -> bool:
        # The analyses left to their default by a study that never reached
        # its analysis: null there is not a choice of none, and was listed
        # against the default list the other study recorded.
        return d.setting == "analysis.include" and any(
            value is None and not side["analysed"]
            for value, side in ((d.first, one), (d.second, other)))

    def as_ran(d: Any) -> Any:
        # Left to the default by a study that ran its analyses: what it ran
        # is what that default was. Listed as "not set" against the same
        # list written out, or against another list.
        if d.setting != "analysis.include":
            return d
        first = one["ran_analyses"] if d.first is None else d.first
        second = other["ran_analyses"] if d.second is None else d.second
        if (isinstance(first, list) and isinstance(second, list)
                and sorted(map(str, first)) == sorted(map(str, second))):
            return None
        return type(d)(d.setting, d.first if first is None else first,
                       d.second if second is None else second)

    found = [d for d in (as_ran(d) for d in found) if d is not None]
    settings = [d.as_record() for d in found if not (d.unrecorded or not_reached(d))]
    unrecorded = [d.setting for d in found if d.unrecorded or not_reached(d)]

    mine, theirs = one["means"], other["means"]
    # An analysis whose record was not read is said so on its side, never
    # as one that recorded no mean.
    for side, means in ((one, mine), (other, theirs)):
        for name in side.get("not_read") or ():
            means.setdefault(name, {"label": _label_of(name), "mean": None, "error": None,
                                    "unit": "", "withheld": None, "not_read": True})
    from fastmdxplora.gui.overview_view import THERMODYNAMICS

    thermodynamic = [key for key, _, _ in THERMODYNAMICS]
    names = set(mine) | set(theirs)
    measures = []
    for name in (_in_order(names - set(thermodynamic))
                 + [key for key in thermodynamic if key in names]):
        first, second = mine.get(name), theirs.get(name)
        row: dict[str, Any] = {"analysis": name, "first": first, "second": second,
                               "label": (first or second or {}).get("label", name),
                               "unit": (first or second or {}).get("unit", "")}
        usable = [side for side in (first, second)
                  if side and side.get("mean") is not None and side.get("error") is not None]
        if len(usable) == 2:
            difference = second["mean"] - first["mean"]
            given = (errors or {}).get(name)
            error = math.hypot(*(given or (first["error"], second["error"])))
            row["versus"] = {"difference": difference, "error": error,
                             "resolved": error > 0 and abs(difference) > RESOLVED_AT * error}
            if given:
                row["versus"]["errors"] = list(given)
                row["versus"]["error_of"] = "replicas"
        measures.append(row)
    return {"settings": settings, "unrecorded": unrecorded, "measures": measures,
            "resolved_at": RESOLVED_AT}


def _phases_ran(base: Path) -> set[str] | None:
    """The phases the study's records say it ran: its Manifest's, and setup
    and simulation where its live record says it simulated. None where
    neither says (its Config's phases are read then)."""
    manifest = _read_json(base / "manifest.json")
    phases = manifest.get("phases") if isinstance(manifest, dict) else None
    ran = {str(p.get("name")) for p in phases or [] if isinstance(p, dict) and p.get("name")
           and str(p.get("status") or "").lower() != "skipped"}
    if _simulated(base):
        ran |= {"setup", "simulation"}
    return ran or None


def _ran_analyses(base: Path) -> list[str] | None:
    """The analyses the study ran, as its analysis manifest planned them."""
    manifest = _read_json(base / "analysis" / "analysis_manifest.json")
    plan = manifest.get("plan") if isinstance(manifest, dict) else None
    return [str(name) for name in plan] if isinstance(plan, list) and _analysed(base) else None


def _analysed(base: Path) -> bool:
    """Whether the study's analysis ran: its Manifest names the phase."""
    manifest = _read_json(base / "manifest.json")
    phases = manifest.get("phases") if isinstance(manifest, dict) else None
    return any(isinstance(p, dict) and p.get("name") == "analysis" for p in phases or [])


# ---------------------------------------------------------------------------
def _means(base: Path) -> list[tuple[str, dict[str, Any]]]:
    """(analysis, {mean, error, unit, withheld}) for each mean the study
    recorded, the measures a card names first leading."""
    from fastmdxplora.batch.aggregate import read_member_findings
    from fastmdxplora.gui.report_dashboard import unit_of

    try:
        found = read_member_findings(base)
    except Exception:  # noqa: BLE001 - a card stands without them
        return []
    out = []
    for name in _in_order(found):
        record = found.get(name) or {}
        mean = record.get("mean") if isinstance(record.get("mean"), dict) else None
        if mean is None:
            continue
        withheld = mean.get("not_a_measurement") or None
        out.append((name, {
            "label": _label_of(name),
            "mean": _finite(mean.get("mean")),
            "error": None if withheld else _finite(mean.get("standard_error")),
            "unit": unit_of(name, mean),
            "withheld": str(withheld) if withheld else None,
        }))
    return out


def _label_of(name: str) -> str:
    """An analysis by the name its series has, else its Analysis page
    heading: "end to end" was its folder's name."""
    from fastmdxplora.gui.report_dashboard import ANALYSIS_SECTION_BY_FOLDER
    from fastmdxplora.gui.series import SERIES

    if name in SERIES:
        return SERIES[name][0]
    return ANALYSIS_SECTION_BY_FOLDER.get(name) or name.replace("_", " ").capitalize()


def _in_order(names: Any) -> list[str]:
    names = list(names)
    return sorted(names, key=lambda n: (_FIRST.index(n) if n in _FIRST else len(_FIRST), n))


def _system_of(config: Any, batch: Any, manifest: Any) -> str:
    for source in (config, batch):
        systems = source.get("systems") if isinstance(source, dict) else None
        if isinstance(systems, list) and systems and isinstance(systems[0], dict):
            named = {str(s.get("system")) for s in systems
                     if isinstance(s, dict) and s.get("system")}
            if len(named) == 1:
                return _short(next(iter(named)))
            if named:
                return f"{len(named)} systems"
    if isinstance(manifest, dict) and manifest.get("system"):
        return _short(str(manifest["system"]))
    return ""


def _short(system: str) -> str:
    """A structure by its identifier, or a file by its name."""
    return Path(system).name if ("/" in system or "\\" in system) else system


def _simulated(base: Path) -> bool:
    """Whether the study's own live record says it simulated: its
    production taken, or steps planned. Analyze again wrote its config's
    phases as `[analysis, report]`, and a run read as a trajectory analysed."""
    record = _read_json(base / "simulation" / "live_status.json")
    if not isinstance(record, dict):
        return False
    states = record.get("stage_states") if isinstance(record.get("stage_states"), dict) else {}
    taken = str(states.get("production") or "waiting").lower() not in ("waiting", "skipped")
    planned = record.get("total_planned_steps")
    return taken or (isinstance(planned, (int, float)) and planned > 0)


def _kind_of(config: Any, batch: Any, *, simulated: bool = False) -> str:
    if isinstance(batch, dict):
        planned = batch.get("planned") or []
        runs = int(batch.get("n_runs") or len(planned) or 0)
        if not runs:
            return "several runs"
        if any(isinstance(p, dict) and ((p.get("options") or {}).get("simulation") or {})
               .get("umbrella") for p in planned):
            return f"umbrella sampling, {runs} windows"
        sweep = batch.get("sweep") if isinstance(batch.get("sweep"), dict) else {}
        axes = sorted(sweep)
        from fastmdxplora.batch.aggregate import SEED_AXES

        if axes and set(axes) <= SEED_AXES:
            return f"{runs} replicas"
        if axes:
            return f"{runs} runs across {', '.join(a.split('.')[-1] for a in axes)}"
        return f"{runs} runs"
    phases = (config or {}).get("include_phase") if isinstance(config, dict) else None
    if isinstance(phases, list) and phases and "simulation" not in phases and not simulated:
        return "a trajectory analysed" if "analysis" in phases else "prepared"
    simulation = (config or {}).get("simulation") if isinstance(config, dict) else None
    if isinstance(simulation, dict):
        for key, said in (("metadynamics", "metadynamics"), ("steered", "steered pulling")):
            if simulation.get(key):
                return said
    return "one run"


def study_state(folder: Path | str) -> str:
    """A study's state, from its records: the one derivation every page
    reads (All studies' cards, the sidebar, the Overview). For a study of
    several runs, from its runs' states (`exploration.runs_of_a_study`)."""
    base = Path(folder)
    return _state_of(base, _read_json(base / "batch_manifest.json"),
                     _read_json(base / "manifest.json"))


def _state_of(base: Path, batch: Any, manifest: Any) -> str:
    try:
        from fastmdxplora.simulation.resume import _still_running

        if _still_running(base):
            return "running"
    except Exception:  # noqa: BLE001 - a record, not a verdict
        pass
    if not isinstance(batch, dict):
        # The live record, before the manifest: the manifest is written when
        # a run ends, so a run cut short has none, or the one before it. A
        # record still saying the run goes on is believed unless the run has
        # ended without saying so (telemetry.status_as_it_stands).
        from fastmdxplora.gui.telemetry import status_as_it_stands

        status = status_as_it_stands(base)
        if status.get("between_rounds"):
            # Its study extends it again in its next round: the Manifest's
            # phases, all ok, say how its pieces so far ended.
            return "running"
        if status.get("recorded_status") and status.get("status") == "completed":
            return "completed"
        said = str(status.get("status") or "").lower()
        if said == "interrupted":
            return "interrupted"
        if said in ("running", "starting", "paused"):
            return "running"
        # A run whose record says it failed or was stopped, as its Overview
        # says: a failure found after the run's own end (the GUI's check of
        # what a run left) is recorded there, and the study was listed
        # Completed beside an Overview saying "Simulation failed".
        if said in ("failed", "stopped"):
            return said
    if isinstance(batch, dict):
        from fastmdxplora.gui.exploration import runs_of_a_study

        runs = runs_of_a_study(base) or []
        states = [run["state"] for run in runs]
        if states and all(state == "completed" for state in states):
            return "completed"
        if any(state == "failed" for state in states):
            return "stopped" if _stopped_in(batch) else "failed"
        # A run stopped or killed in a later piece (its study's rounds).
        if any(state == "stopped" for state in states):
            return "stopped"
        # A run going on, as its record and its process say
        # (`runs_of_a_study` calls one whose process is gone interrupted).
        if any(state == "running" for state in states):
            return "running"
        return "incomplete"
    phases = manifest.get("phases") if isinstance(manifest, dict) else None
    if not isinstance(phases, list) or not phases:
        return "not started"
    for phase in phases:
        if isinstance(phase, dict) and phase.get("status") == "error":
            code = str((phase.get("refusal") or {}).get("code") or "")
            return "stopped" if code == "simulation.run.stopped" else "failed"
    return "completed"


def _stopped_in(batch: dict[str, Any]) -> bool:
    return any(isinstance(run, dict) and run.get("error_type") == "Stopped"
               for run in batch.get("runs") or [])


def _when(base: Path, manifest: Any) -> str:
    phases = manifest.get("phases") if isinstance(manifest, dict) else None
    if isinstance(phases, list):
        for phase in phases:
            if isinstance(phase, dict) and phase.get("started_at"):
                return str(phase["started_at"])
    try:
        stamp = os.path.getmtime(base)
    except OSError:
        return ""
    return datetime.fromtimestamp(stamp, timezone.utc).isoformat(timespec="seconds")


#: The records a study's activity changes: its run, its analyses and report,
#: and what the person keeps with it (tags, views, selections). The Agent's
#: conversations are looked at as a folder, each file in it.
_ACTIVE_RECORDS = (
    "manifest.json", "batch_manifest.json", "simulation/live_status.json",
    "analysis/analysis_manifest.json", "report/report.md", "study_tags.json",
    "viewer_views.json", "viewer_selections.json",
)
_CONVERSATIONS = Path("agent") / "conversations"


def _last_active(base: Path) -> str:
    """When the study was last active, by its newest record, or ''."""
    newest = 0.0
    for name in _ACTIVE_RECORDS:
        try:
            newest = max(newest, (base / name).stat().st_mtime)
        except OSError:
            continue
    try:
        with os.scandir(base / _CONVERSATIONS) as kept:
            for entry in kept:
                if entry.is_file():
                    newest = max(newest, entry.stat().st_mtime)
    except OSError:
        pass
    if not newest:
        return ""
    return datetime.fromtimestamp(newest, timezone.utc).isoformat(timespec="seconds")


def _lifted(config: Any) -> Any:
    """A study of several runs writes its phase settings under `options`;
    lifted beside the rest, so one reading serves both."""
    if isinstance(config, dict) and isinstance(config.get("options"), dict):
        return {**config, **config["options"]}
    return config


def _read_json(path: Path) -> Any:
    # Through `record_access`: held to its study where the read is (a study
    # of runs' state said on `/api/status`), as pathlib elsewhere.
    from fastmdxplora import record_access

    try:
        return json.loads(record_access.read_text(path, encoding="utf-8"))
    except (OSError, ValueError):
        return None


def _read_yaml(path: Path) -> Any:
    import yaml

    from fastmdxplora import record_access

    try:
        return yaml.safe_load(record_access.read_text(path, encoding="utf-8"))
    except (OSError, yaml.YAMLError):
        return None


def _finite(value: Any) -> float | None:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return None
    return float(value) if math.isfinite(value) else None
