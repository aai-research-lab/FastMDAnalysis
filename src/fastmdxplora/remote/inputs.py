"""The files a config points at, gathered so they can travel with it.

A config names files on this computer: the structure, a ligand's SDF, force
field XMLs, a trajectory to analyse, a prepared study to continue from.
None of them exist on the machine the study is sent to. So every string in
the config that names an existing file or directory -- relative to the
config's own folder, as the config means it -- is sent under ``inputs/``,
and the copy of the config that travels names it there.

Read from the values rather than from a list of path-valued settings, so a
setting added tomorrow that takes a file travels without anyone remembering
to add it here.

A structure named by PDB ID is fetched by the machine itself, which needs
the internet; the caller refuses where it has none.

**Only what is in the study's own folder travels.** A file the config names
outside the folder holding it is refused, not sent, whoever wrote the config:
a config naming ``~/.ssh/id_ed25519`` would otherwise copy the key to the
machine. Links are resolved first, so a link in the folder
pointing out of it is outside too, and a folder that travels is refused if
a link in it, or in a folder one of its links leads to, points out of it or
back into itself, since the copy follows links. One
exception: a prepared study named by ``simulation.setup_from`` or
``simulation.resume_from`` may sit beside the folder, where a study usually
does (anywhere in the folder holding the study's folder), when it holds a
manifest FastMDXplora wrote; links in it are held to that study's folder.

Nothing in a place keys and credentials are kept (``.ssh``, ``.gnupg``,
``.aws`` and the like, and FastMDXplora's own settings) travels: a folder
that travels is refused if one is anywhere in it.
"""

from __future__ import annotations

import hashlib
import json
import os
import re
import stat
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from fastmdxplora.refusals import StudyError

__all__ = ["Inputs", "fingerprint_of", "gather_inputs", "link_out_of", "outside_said",
           "kept_place_of", "kept_place_said", "private_in", "refuse_home_or_above",
           "refuse_kept_place", "refuse_no_study", "names_a_study", "size_of"]

#: Settings that name where results go rather than something to read.
_NOT_INPUTS = frozenset({"output"})

#: A four-character PDB identifier, which setup fetches from RCSB.
_PDB_ID = re.compile(r"^[0-9][A-Za-z0-9]{3}$")

#: Settings naming a prepared study, which may sit beside the folder.
_STUDIES_NAMED = frozenset({"setup_from", "prepared_from", "resume_from"})


@dataclass
class Inputs:
    """The config as it travels, and what travels with it."""

    config: Any
    #: name under ``inputs/`` -> the file or directory on this computer.
    files: dict[str, Path] = field(default_factory=dict)
    #: Structures the machine will fetch from RCSB itself.
    fetched: list[str] = field(default_factory=list)
    #: name under ``inputs/`` -> the folder its links must stay inside.
    held_to: dict[str, Path] = field(default_factory=dict)


#: The suffixes setup reads a ``system`` by as a structure file
#: (``setup.pipeline``).
STRUCTURE_SUFFIXES = frozenset({".pdb", ".cif", ".pdbx"})


#: The settings whose values are words, names or selections and never a
#: file, though the schema gives them no list of choices: chains, a
#: mutation, a residue or a water model's name, a selection, a title, a
#: decision's reasons, a system's name.
_WORDS_WITHOUT_CHOICES = frozenset({
    "include_phase", "exclude_phase", "agent_model", "mutations", "mutation_chain",
    "chains", "water_model", "ligand_name", "ligand_resname", "ligand_pose",
    "ion_positive", "ion_negative", "constraints", "device_index", "devices",
    "save_selection", "dashboard_ligand_resname", "restrain", "select_atoms",
    "selection", "title", "author", "region_highlights", "id", "why", "source",
    "alternatives", "selection_a", "selection_b", "select_atoms_a", "select_atoms_b",
    "bilayer_selection", "site_selection", "axis_selection"})


def word_keys() -> frozenset[str]:
    """The settings whose values are words the software knows, never a file:
    each with a list of choices in the schema (the analyses, a force field's
    name, a box's shape), the phases a study runs, and the names and
    selections of ``_WORDS_WITHOUT_CHOICES``."""
    from fastmdxplora.config.schema import all_schemas

    return _WORDS_WITHOUT_CHOICES | frozenset(
        field.name for schema in all_schemas().values() for field in schema.fields
        if field.choices)


def read_as_a_name(system: str) -> bool:
    """Whether setup reads ``system`` as a PDB ID or a sequence, as
    ``setup.pipeline`` classifies it (no structure suffix; four letters and
    digits, or letters alone), wherever it runs: never a file."""
    try:
        if Path(system).suffix.lower() in STRUCTURE_SUFFIXES:
            return False
    except (TypeError, ValueError):
        return False
    return (len(system) == 4 and system.isalnum()) or system.isalpha()


def _existing(value: str, base: Path) -> Path | None:
    if not value or "\n" in value or len(value) > 4096:
        return None
    try:
        # "~100 ns" is words, not a home folder: a home that cannot be found
        # makes it no path.
        path = Path(value).expanduser()
    except RuntimeError:
        return None
    if not path.is_absolute():
        path = base / path
    try:
        return path.resolve() if path.exists() else None
    except OSError:
        return None


def _inside(path: Path, folder: Path) -> bool:
    return path == folder or folder in path.parents


def _walked(path: Path):
    """Each folder and file under ``path`` as a copy that follows links
    reaches it: the place it is reached at, and where it really is (``None``
    where a link leads nowhere, or back into a folder it is inside, which
    the copy would follow without end)."""
    within: dict[Path, tuple[Path, ...]] = {}
    unread: list[str] = []
    # A folder that cannot be listed cannot be checked: said as unresolved.
    for top, dirs, files in os.walk(path, followlinks=True,
                                    onerror=lambda e: unread.append(e.filename)):
        while unread:
            yield Path(unread.pop()), None
        here = Path(top)
        real = here.resolve()
        above = within.get(here.parent, ()) if here != path else ()
        if real in above:
            dirs[:] = []
            yield here, None
            continue
        within[here] = (*above, real)
        yield here, real
        for name in files:
            entry = here / name
            try:
                yield entry, entry.resolve(strict=True)
            except (OSError, RuntimeError):
                yield entry, None
        for name in list(dirs):
            entry = here / name
            if entry.is_symlink() and not entry.exists():
                dirs.remove(name)
                yield entry, None
    while unread:
        yield Path(unread.pop()), None


def link_out_of(path: Path, folder: Path) -> Path | None:
    """The first place in ``path`` a copy that follows links would be led
    out of ``folder`` by, if any.

    ``path`` is a file or a folder. Its links are followed as the copy
    follows them, through a folder a link leads to and the links in it in
    turn, and every file and folder reached is resolved: one outside
    ``folder``, one a link leads nowhere, or a folder reached twice (a
    link back into itself, which the copy would follow without end) is
    returned.
    """
    try:
        real = path.resolve(strict=True)
    except (OSError, RuntimeError):
        return path
    if not _inside(real, folder):
        return path
    if not path.is_dir():
        return None
    for entry, real in _walked(path):
        if real is None or not _inside(real, folder):
            return entry
    return None


#: Where keys and credentials are kept: never sent, whatever folder holds them.
_PRIVATE_PARTS = frozenset({".ssh", ".gnupg", ".aws", ".netrc", ".kube", ".docker",
                            ".pgpass", ".git-credentials", ".azure", ".password-store",
                            ".vault-token", ".pypirc", ".npmrc"})
#: The same, as a folder inside another.
_PRIVATE_PAIRS = frozenset({(".config", "gcloud"), (".config", "gh"),
                            (".config", "rclone"), ("library", "keychains")})


#: FastMDXplora's own settings: the Agent's key, the machines, the jobs.
_SETTINGS_KEPT = ("model.json", "calibration.json", "machines", "jobs", "sockets",
                  "gpu_memory", "gpu_memory_here.json")


def _settings_kept() -> tuple[str, ...]:
    """Where settings and keys are kept, each resolved: never sent.
    FastMDXplora's own (:func:`own_settings_kept`; its settings folder
    itself may hold a study, only these may not), and other programs' by
    name or by a link (:func:`_kept_by_name_or_link`)."""
    return own_settings_kept() + _kept_by_name_or_link()


def _own_settings(places: list[Path]) -> tuple[str, ...]:
    """FastMDXplora's own settings in each of its settings folders."""
    return tuple(str(settings / name) for settings in dict.fromkeys(places)
                 for name in _SETTINGS_KEPT)


def own_settings_kept() -> tuple[str, ...]:
    """FastMDXplora's own settings alone (its key, machines, jobs): the
    places of keys :func:`_settings_kept` holds, without other programs'."""
    from fastmdxplora.user_dir import user_config_dir

    places = []
    try:
        places.append(user_config_dir().resolve())
    except (OSError, RuntimeError, ValueError):
        pass
    for home in home_folders():
        try:
            places.append((home / ".config" / "fastmdxplora").resolve())
        except (OSError, RuntimeError, ValueError):
            continue
    return _own_settings(places)


#: Folders in the home folder's ``Library`` that hold the person's own files.
_LIBRARY_FILES = frozenset({"cloudstorage", "mobile documents"})
#: How deep links are looked for in each place keys or settings are kept,
#: and how many entries one place may hold before the search gives up and
#: nothing is sent.
_LINK_DEPTH = 3
_LINK_ENTRIES = 20_000

#: Folders not looked in for links: what was thrown away.
_NOT_LOOKED_IN = frozenset({"trash", ".trash"})

_KEPT_CACHE: dict[tuple, tuple[float, tuple[str, ...]]] = {}
#: Where each place found by a link was found from, to say so.
_LINKED_FROM: dict[str, str] = {}


def _places_looked_in(homes: tuple[Path, ...], named: tuple[str, ...]
                      ) -> tuple[list[Path], list[Path]]:
    """The places settings are kept by name that are not named by a dot
    (the folders the XDG and APPDATA variables name, where full paths, and
    each home's ``AppData``), and the places links are looked for in: each
    home's places of keys (:data:`_PRIVATE_PARTS`, ``.config`` with its
    pairs, ``Library/Keychains``), its ``.config``, ``.local/share`` and
    those named folders. Not ``Library/Containers``, Wine's folders, the
    Trash or caches: their links lead to the person's own folders."""
    by_name: list[Path] = []
    for text in named:
        try:
            place = Path(os.path.normpath(Path(text).expanduser()))
        except (RuntimeError, ValueError):
            continue
        if text and place.is_absolute():
            by_name.append(place)
    looked: list[Path] = list(by_name)
    for home in homes:
        by_name.append(home / "AppData")
        looked += [home / name for name in sorted(_PRIVATE_PARTS)]
        looked += [home / ".config", home / ".local" / "share",
                   home / "Library" / "Keychains"]
    return by_name, looked


def _kept_by_name_or_link() -> tuple[str, ...]:
    """Where settings and keys are kept that no dot directly in a home names:
    the folders of :func:`_places_looked_in` kept by name, and wherever a
    link leads from a home's own entry whose name starts with a dot (``~/.ssh
    -> ~/dotfiles/ssh``) or from up to :data:`_LINK_DEPTH` folders inside a
    place looked in (``~/.config/gh -> ~/dotfiles/gh``, a key a dotfiles
    manager links into ``~/.ssh``), each resolved. A link to a home folder,
    a folder above one, or one of the person's own folders directly in a
    home (``~/Documents``) is no place of keys. A folder that cannot be read
    is passed over; a place with more than :data:`_LINK_ENTRIES` entries is
    refused, as where its links lead cannot be told. Looked up again after
    ten seconds at most."""
    homes = tuple(home_folders())
    named = tuple(os.environ.get(name, "") for name in
                  ("XDG_CONFIG_HOME", "XDG_DATA_HOME", "APPDATA", "LOCALAPPDATA"))
    key = (homes, named)
    now = time.monotonic()
    kept = _KEPT_CACHE.get(key)
    if kept is not None and now - kept[0] < 10:
        return kept[1]
    above = {Path(p.anchor) for p in homes}
    for home in homes:
        above.update((home, *home.parents))
    by_name, looked = _places_looked_in(homes, named)
    links: list[Path] = []
    for home in homes:
        try:
            with os.scandir(home) as listed:
                links += [Path(entry.path) for entry in listed
                          if entry.name.startswith(".") and entry.is_symlink()]
        except OSError:
            pass
    for place in looked:
        left = _LINK_ENTRIES
        stack = [(place, 0)]
        while stack:
            folder, depth = stack.pop()
            try:
                with os.scandir(folder) as listed:
                    for entry in listed:
                        left -= 1
                        if left < 0:
                            raise StudyError(
                                f"{place} holds more than {_LINK_ENTRIES} entries, so "
                                "where links in it lead cannot be told, and nothing is "
                                "sent. Move what is not settings out of it.",
                                code="remote.input.outside", given=str(place),
                                where="", folder=str(place))
                        if entry.is_symlink():
                            links.append(Path(entry.path))
                        elif (depth + 1 < _LINK_DEPTH and entry.name.casefold()
                              not in _NOT_LOOKED_IN and entry.is_dir(follow_symlinks=False)):
                            stack.append((Path(entry.path), depth + 1))
            except OSError:
                continue
    places: list[str] = []
    for place in by_name:
        if place not in above:
            places.append(str(place))
    _LINKED_FROM.clear()
    for link in links:
        try:
            target = link.resolve()
        except (OSError, RuntimeError, ValueError):
            continue
        # The person's own folders are not made places of keys by a link.
        if target in above or (target.parent in homes and not target.name.startswith(".")):
            continue
        places.append(str(target))
        _LINKED_FROM.setdefault(str(target), str(link))
    kept = tuple(dict.fromkeys(places))
    _KEPT_CACHE.clear()
    _KEPT_CACHE[key] = (now, kept)
    return kept


def _private(path: Path, kept: tuple[str, ...] = ()) -> bool:
    parts = [part.casefold() for part in path.parts]
    if _PRIVATE_PARTS & set(parts) or any(
            pair in _PRIVATE_PAIRS for pair in zip(parts, parts[1:])):
        return True
    # By the path's text: a folder's parents are slow to walk 20,000 times.
    # Where the disk ignores the case of letters, so does this. Each folder
    # above it looked up among the places, not each place tried in turn.
    places = _places_as_read(kept)
    text = _as_the_disk_reads(str(path))
    while True:
        if text in places:
            return True
        above = os.path.dirname(text)
        if above == text:
            return False
        text = above


_PLACES_READ: dict[tuple[str, ...], frozenset[str]] = {}


def _places_as_read(kept: tuple[str, ...]) -> frozenset[str]:
    """``kept`` as the disk compares names, as a set (the last kept)."""
    found = _PLACES_READ.get(kept)
    if found is None:
        _PLACES_READ.clear()
        found = _PLACES_READ[kept] = frozenset(
            os.path.normpath(_as_the_disk_reads(place)) for place in kept)
    return found


def _as_the_disk_reads(text: str) -> str:
    """``text`` as this platform's usual disk compares names: macOS and
    Windows ignore the case of letters."""
    return text.casefold() if sys.platform in ("darwin", "win32") else text


def private_in(path: Path) -> Path | None:
    """The first place in ``path`` keys and credentials are kept, if any:
    ``path`` itself, or a file or folder a copy that follows links reaches
    in it, by where it is reached and where it really is."""
    kept = _settings_kept()
    try:
        real = path.resolve()
    except (OSError, RuntimeError):
        real = path
    if _private(path, kept) or _private(real, kept):
        return path
    if not path.is_dir():
        return None
    for entry, real in _walked(path):
        if _private(entry, kept) or (real is not None and _private(real, kept)):
            return entry
    return None


def _home_or_above(folder: Path) -> bool:
    """Whether ``folder`` is the top of the file system, the home folder or
    a folder holding it, also where the disk ignores the case of letters;
    True where no home folder can be found, since then it cannot be told."""
    return home_or_above(folder, home_unknown=True)


def home_folders() -> list[Path]:
    """The home folder as the account's own entry names it, then as
    ``HOME`` names it, each resolved, or none where neither can be found:
    a wrong or empty ``HOME`` does not hide the account's home, and a
    ``HOME`` that is not a full path, or is a disk's top, names none."""
    found: list[Path] = []
    try:
        import pwd

        entry = pwd.getpwuid(os.getuid()).pw_dir
        if entry and Path(entry).is_absolute():
            found.append(Path(entry).resolve())
    except (ImportError, AttributeError, KeyError, OSError, RuntimeError, ValueError):
        pass
    try:
        named = Path.home()
        if named.is_absolute():
            found.append(named.resolve())
    except (RuntimeError, OSError, KeyError, ValueError):
        pass
    return [home for home in dict.fromkeys(found) if home != Path(home.anchor)]


def home_or_above(folder: Path, *, home_unknown: bool = False) -> bool:
    """As :func:`_home_or_above`; ``home_unknown`` is the answer where the
    home folder cannot be found (no ``HOME`` and no account entry)."""
    if folder == Path(folder.anchor):
        return True
    homes = home_folders()
    if not homes:
        return home_unknown
    for home in homes:
        for candidate in (home, *home.parents):
            if folder == candidate:
                return True
            try:
                if os.path.samefile(folder, candidate):
                    return True
            except OSError:
                continue
    return False


def size_of(path: Path) -> int:
    """Bytes in a file, or in the files a copy that follows links sends of
    a folder."""
    try:
        if not path.is_dir():
            return path.stat().st_size
        return sum(real.stat().st_size for _, real in _walked(path)
                   if real is not None and real.is_file())
    except OSError:
        return 0


#: Files up to this size are read whole for a send's fingerprint, until
#: :data:`READ_IN_ALL_BYTES` have been read for one send; past either, a
#: file is known by its size, the file it is (device and inode) and when it
#: was last written and changed (``ctime``, which a copy that keeps times
#: does not put back).
READ_WHOLE_BYTES = 256 * 1024 * 1024
READ_IN_ALL_BYTES = 2 * 1024 * 1024 * 1024


def _file_print(real: Path, left: list[int]) -> tuple[int, str]:
    """A file's size, and what tells it from the same file changed;
    ``left`` is what may still be read for this send, and is spent."""
    # Opened without waiting, and what was opened is what is judged: a file
    # put in place of a pipe between a look and an open never holds the
    # reader, and a pipe or a device, which a copy does not send, is never
    # read.
    handle = os.open(real, os.O_RDONLY | getattr(os, "O_NONBLOCK", 0)
                     | getattr(os, "O_BINARY", 0))
    with os.fdopen(handle, "rb") as opened:
        found = os.fstat(opened.fileno())
        if not stat.S_ISREG(found.st_mode):
            return 0, f"special:{stat.S_IFMT(found.st_mode)}"
        if found.st_size > READ_WHOLE_BYTES or found.st_size > left[0]:
            return found.st_size, (f"{found.st_size}:{found.st_dev}:{found.st_ino}:"
                                   f"{found.st_mtime_ns}:{found.st_ctime_ns}")
        left[0] -= found.st_size
        digest = hashlib.sha256()
        while chunk := opened.read(1 << 20):
            digest.update(chunk)
    return found.st_size, f"{found.st_size}:{digest.hexdigest()}"


def fingerprint_of(path: Path, left: list[int] | None = None) -> tuple[int, str]:
    """The bytes a copy that follows links sends of ``path`` (as
    :func:`size_of`), and a digest that changes when anything it sends
    does: a file's contents (or, past what is read whole, the file it is and
    when it was written and changed), and in a folder each file's place in
    it. ``left`` is what may still be read for one send, shared by its
    inputs (:data:`READ_IN_ALL_BYTES` where not given)."""
    left = [READ_IN_ALL_BYTES] if left is None else left
    digest = hashlib.sha256()
    total = 0
    if not path.is_dir():
        try:
            total, said = _file_print(path.resolve(strict=True), left)
        except (OSError, RuntimeError, ValueError):
            total, said = 0, "unread"
        digest.update(said.encode())
        return total, digest.hexdigest()
    lines = []
    try:
        # In one order whatever order the folder lists them in, so the
        # files read whole, before what may be read is spent, are the same
        # at the plan and at the send.
        walked = sorted(_walked(path), key=lambda entry: entry[0].as_posix())
    except (OSError, RuntimeError, ValueError):
        walked = [(path, None)]
    for here, real in walked:
        place = here.relative_to(path).as_posix() if here.is_relative_to(path) else str(here)
        said = "unread"
        if real is not None:
            try:
                if real.is_dir():
                    said = "folder"
                else:
                    size, said = _file_print(real, left)
                    total += size
            except (OSError, RuntimeError, ValueError):
                said = "unread"
        lines.append(f"{place}\0{said}\n")
    digest.update("".join(sorted(lines)).encode())
    return total, digest.hexdigest()


def _study_folder(path: Path) -> Path | None:
    """The study ``path`` is or is the setup folder of, by its manifest:
    one FastMDXplora wrote, which lists the phases it ran."""
    for folder in (path, path.parent) if path.name == "setup" else (path,):
        try:
            record = json.loads((folder / "manifest.json").read_text(encoding="utf-8"))
        except (OSError, ValueError):
            continue
        if isinstance(record, dict) and isinstance(record.get("phases"), list):
            return folder
    return None


def outside_said(where: str, given: str, found: Path, folder: Path,
                 *, link: Path | None = None) -> str:
    """Why a file named outside the study's folder is not sent."""
    if link is not None:
        said = (f"`{where}` names {given}, and {link} in it is a link "
                f"leading out of {folder}, so it is not sent: what it points "
                "at would travel with the study. Copy the file in its place.")
    else:
        said = (f"`{where}` names {given}, which is {found}, outside the "
                f"study's folder {folder}, so it is not sent. Only files in "
                "the folder holding the config travel with it. Copy it into "
                "that folder and name it there"
                + (", or name a structure by its PDB identifier"
                   if where.rsplit(".", 1)[-1].split("[")[0] in ("system", "systems")
                   else "")
                + ".")
    return said


def kept_place_of(path: Path, *, a_file: bool = False) -> Path | None:
    """The place keys or settings are kept that ``path`` is in, as named or
    where its links lead, or None: one of :data:`_PRIVATE_PARTS` or the
    pairs, FastMDXplora's own settings, a folder settings are kept in by
    name or by a link (:func:`_kept_by_name_or_link`), or a folder directly
    in a home folder whose name starts with a dot (``~/.config``,
    ``~/.local``), or that is ``AppData`` or ``Library`` (its
    ``CloudStorage`` and ``Mobile Documents`` aside); of ``a_file``, its
    folders only. A study's config there, or its results, never travel."""
    kept = _settings_kept()
    places = [Path(os.path.normpath(Path(path).expanduser()))]
    try:
        places.append(Path(path).expanduser().resolve())
    except (OSError, RuntimeError, ValueError):
        pass
    homes = home_folders()
    for place in places:
        # The file itself too: a single file a link in a place of keys leads to.
        parts = place.parts
        for depth in range(1, len(parts) + 1):
            if _private(Path(*parts[:depth]), kept):
                return Path(*parts[:depth])
        folded = [_as_the_disk_reads(part) for part in place.parts]
        for home in homes:
            ours = [_as_the_disk_reads(part) for part in home.parts]
            if folded[:len(ours)] != ours or len(folded) == len(ours):
                continue
            rest = place.parts[len(ours):]
            if len(rest) < (2 if a_file else 1):
                continue
            first = rest[0]
            if first.startswith(".") or first.casefold() == "appdata" or (
                    first.casefold() == "library"
                    and (len(rest) < 2 or rest[1].casefold() not in _LIBRARY_FILES)):
                return home / first
    return None


def refuse_kept_place(path: Path, *, results: bool = False,
                      output_given: bool = False) -> None:
    """Refuse (``remote.input.outside``) a config kept in a place keys or
    settings are kept (:func:`kept_place_of`), which would travel as the
    study's config, or a results folder there, which a fetch would fill;
    ``output_given`` where `--output` named it."""
    place = kept_place_of(Path(path), a_file=not results)
    if place is None:
        return
    raise StudyError(kept_place_said(path, place, results=results,
                                     output_given=output_given),
                     code="remote.input.outside", given=str(path), where="",
                     folder=str(place))


def names_a_study(raw: Any) -> bool:
    """Whether a config's text names a study: its `systems` (or `system`),
    or the study it carries on (`simulation.resume_from` or `setup_from`)."""
    if not isinstance(raw, dict):
        return False
    if raw.get("systems") or raw.get("system"):
        return True
    simulation = raw.get("simulation")
    return isinstance(simulation, dict) and bool(
        simulation.get("resume_from") or simulation.get("setup_from"))


def refuse_no_study(raw: Any, config_path: Path) -> None:
    """Refuse a file that names no study (:func:`names_a_study`): no
    study's config, so a file of settings or keys never travels as one,
    wherever it is kept."""
    if names_a_study(raw):
        return
    raise StudyError(
        f"{config_path} names no `systems`, nor a study it carries on, so it is no "
        "study's config and it is not sent. Give the config that names the study's "
        "systems.",
        code="config.option.missing_companion", option="systems", requires=["systems"])


def kept_place_said(path: Any, place: Path, *, results: bool = False,
                    in_gui: bool = False, output_given: bool = False,
                    a_study: bool = True) -> str:
    """Why a config (``a_study``; else a file) or a results folder in a place
    keys or settings are kept is refused, with what to do instead, where it
    is said: at the command line `--output` names the results folder (or,
    not given, the folder the send is run from); in the GUI, `output` in
    the config or **Results folder** in the form does. Where the place was
    found by a link, the link is named."""
    found = _LINKED_FROM.get(str(place))
    where = f"{place} (a link, {found}, leads there)" if found else str(place)
    if results:
        if in_gui:
            how = ("Set `output` in the config (or Results folder in the form) to a "
                   "folder outside such places.")
        elif output_given:
            how = "Give `--output` naming a folder outside such places."
        else:
            how = ("Give `--output` naming a folder outside such places, or run the "
                   "send from another folder.")
        return (f"The results folder {path} is in {where}, where keys or settings are "
                f"kept, so nothing is fetched into it. {how}")
    what = "config" if a_study else "file"
    return (f"The {what} {path} is in {where}, where keys or settings are kept, so it "
            "is not sent. Put the study's config and its files in a folder of their "
            "own.")


def refuse_home_or_above(folder: Path) -> None:
    """Refuse (``remote.input.outside``) a config whose folder, which is
    what travels with it, is the top of the file system or the home folder,
    or where no home folder can be found, so that cannot be told."""
    folder = folder.resolve()
    if not _home_or_above(folder):
        return
    said = (f"The config is in {folder}, and a study's folder is what travels with "
            "it; that cannot be the top of the file system or your home folder. "
            "Put the config and its files in a folder of their own.")
    if folder != Path(folder.anchor) and not home_folders():
        said = (f"The config is in {folder}, and a study's folder is what travels "
                "with it. No home folder was found, so FastMDXplora cannot tell whether "
                "that is your home folder, which never travels. Set `HOME` to your "
                "home folder and send it again.")
    raise StudyError(said, code="remote.input.outside", given=str(folder), where="",
                     folder=str(folder))


def gather_inputs(config: Any, base: Path) -> Inputs:
    """``config`` with local paths renamed under ``inputs/``, and the files.

    Refuses (``remote.input.outside``) a file outside ``base`` or a folder
    with a link leading out of it, except a prepared study named by
    ``setup_from`` or ``resume_from``, as the module says.
    """
    folder = base.resolve()
    refuse_home_or_above(folder)
    found = Inputs(config=None)
    by_source: dict[Path, str] = {}

    def travel(path: Path, held_to: Path) -> str:
        if path in by_source:
            return by_source[path]
        name = path.name or "input"
        taken = set(found.files)
        candidate, n = name, 1
        while candidate in taken:
            n += 1
            candidate = f"{n}-{name}"
        found.files[candidate] = path
        found.held_to[candidate] = held_to
        by_source[path] = candidate
        return candidate

    def allowed(path: Path, key: str, where: str, given: str) -> Path:
        """The folder ``path``'s links are held to, or a refusal."""
        held_to = folder if _inside(path, folder) else None
        if (held_to is None and key in _STUDIES_NAMED and path.is_dir()
                and where.endswith(f"simulation.{key}")
                and folder.parent in path.parents):
            held_to = _study_folder(path)
            if held_to is not None and _home_or_above(held_to):
                held_to = None
        link = None if held_to is None else link_out_of(path, held_to)
        kept = None if held_to is None or link is not None else private_in(path)
        if kept is None and held_to is not None and link is None:
            # A prepared study kept where programs keep settings.
            kept = kept_place_of(path)
        if kept is not None:
            raise StudyError(
                f"`{where}` names {given}, "
                + ("which is" if kept == path else f"and {kept} in it is")
                + " a place keys and credentials are kept, so it is not sent.",
                code="remote.input.outside", given=given, where=where,
                folder=str(held_to))
        if held_to is None or link is not None:
            refused_in = held_to or folder
            raise StudyError(
                outside_said(where, given, path, refused_in, link=link),
                code="remote.input.outside", given=given, where=where,
                folder=str(refused_in))
        return held_to

    words = word_keys()

    def walk(value: Any, key: str = "", where: str = "") -> Any:
        if key in words:
            # "setup", "rmsd": a phase or an analysis, though a folder of
            # that name sits beside the config (a study's own folder).
            return value
        if isinstance(value, dict):
            return {k: (v if k in _NOT_INPUTS and not where
                        else walk(v, k, f"{where}.{k}" if where else str(k)))
                    for k, v in value.items()}
        if isinstance(value, list):
            return [walk(item, key, f"{where}[{i}]") for i, item in enumerate(value)]
        if isinstance(value, str):
            if key in ("system", "systems") and read_as_a_name(value):
                # A PDB ID or a sequence to setup, wherever it runs: never a
                # file or folder of that name beside the config.
                if _PDB_ID.match(value):
                    found.fetched.append(value)
                return value
            path = _existing(value, folder)
            if path is not None:
                return f"inputs/{travel(path, allowed(path, key, where, value))}"
        return value

    found.config = walk(config)
    return found
