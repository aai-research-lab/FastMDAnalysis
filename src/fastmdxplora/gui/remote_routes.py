"""The GUI's way to your machines: the routes behind **Run on** and
**Remote jobs**.

What ``fastmdx remote`` does at a terminal, for the page, under the rules
every interface keeps (`docs/remote.md`):

- machines only from the records made at a terminal; nothing here inspects a
  machine or installs on one;
- a config, the files that travel and the results folder all inside the
  GUI's workspace, and only the files in the config's folder travel;
- nothing sent until the person agrees to the plan the page showed: a plan
  is kept here under a token for ten minutes and used once, and the send
  goes ahead only where what would travel now is what was shown, each file's
  contents included (its digest, as an AI app's yes is bound to it);
- a fetch only of the size the person was shown;
- refused in a hosted GUI, whose runs are its service's.

The server's gate holds for every route here as for the others: loopback
only, ``Host`` and ``Origin`` checked.
"""

from __future__ import annotations

import secrets
import threading
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from fastmdxplora.refusals import StudyError, refusal_of

__all__ = ["RemoteDesk", "PLAN_KEPT_S"]

#: How long a plan shown on the page may be sent.
PLAN_KEPT_S = 600

#: Plans kept at once; the oldest goes first.
_MOST_PLANS = 32


@dataclass(frozen=True)
class _Plan:
    config: Path
    machine: str
    output: Path
    digest: str
    made_at: float


def _read_whole(path: Path) -> tuple[str, int, int] | None:
    """The sha256 of a regular file's bytes, with its device and inode, read
    through one handle that follows no link and waits on no pipe; None
    where it is no such file or cannot be read."""
    import hashlib
    import os
    import stat

    flags = (os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0) | getattr(os, "O_NONBLOCK", 0)
             | getattr(os, "O_BINARY", 0))
    try:
        handle = os.open(path, flags)
    except OSError:
        return None
    try:
        found = os.fstat(handle)
        if not stat.S_ISREG(found.st_mode):
            return None
        digest = hashlib.sha256()
        while True:
            block = os.read(handle, 1 << 20)
            if not block:
                break
            digest.update(block)
        return digest.hexdigest(), found.st_dev, found.st_ino
    except OSError:
        return None
    finally:
        os.close(handle)


def _sha_of(path: Path) -> str:
    """The sha256 of a regular file's bytes, or "" where it cannot be read."""
    read = _read_whole(path)
    return read[0] if read else ""


def _taken_back(path: Path, sha: str) -> None:
    """Delete ``path`` only while it is the regular file of that sha256 it
    was written as: never a link, a file changed or put in its place."""
    import os

    read = _read_whole(path)
    if not sha or read is None or read[0] != sha:
        return
    try:
        now = os.lstat(path)
        if (now.st_dev, now.st_ino) == read[1:]:
            path.unlink()
    except OSError:
        pass


def _refused(error: str, code: str = "", **more: Any) -> dict[str, Any]:
    return {"ok": False, "error": error, **({"code": code} if code else {}), **more}


#: The routes, so another path under ``/api/remote/`` is not found.
_GETS = frozenset({"/api/remote/machines", "/api/remote/job", "/api/remote/fetch-sizes"})
_POSTS = frozenset({"/api/remote/plan", "/api/remote/send", "/api/remote/forget",
                    "/api/remote/fetch",
                    "/api/remote/cancel"})

_HOSTED = _refused("Not available in a hosted GUI: its studies run on its service.",
                   hosted=True)


def _read_otherwise(text: str, sent_from: Path, run_from: Path) -> list[str]:
    """Each value in a config's text that names, by a relative path, a file
    or folder a send would carry from ``sent_from`` and a run here read from
    ``run_from`` as another or not at all, or a path (a name with a folder or
    a suffix) a run here would read and a send would not carry: the results
    folder aside, and a ``system`` that setup reads as a PDB ID or a
    sequence, and a setting whose values are words (a phase, an analysis).
    Any other word is left as a word."""
    import yaml

    from fastmdxplora.remote.inputs import _existing, read_as_a_name, word_keys

    try:
        config = yaml.safe_load(text)
    except yaml.YAMLError:
        return []
    found: list[str] = []
    words = word_keys()

    def walk(value: Any, key: str, top: bool) -> None:
        if key in words:
            return
        if isinstance(value, dict):
            for name, item in value.items():
                if not (top and name == "output"):
                    walk(item, str(name), False)
        elif isinstance(value, list):
            for item in value:
                walk(item, key, False)
        elif isinstance(value, str):
            if key in ("system", "systems") and read_as_a_name(value):
                return
            try:
                relative = not Path(value).expanduser().is_absolute()
            except (RuntimeError, ValueError):
                return
            if not relative:
                return
            sent, run = _existing(value, sent_from), _existing(value, run_from)
            # What a send would carry, read otherwise here; or a path that a
            # run here reads and a send would not carry at all.
            if sent != run and (sent is not None or "/" in value or "\\" in value
                                or Path(value).suffix):
                found.append(value)

    walk(config, "", True)
    return found


def _absolute(given: Any) -> bool:
    """Whether a path as named is a full one (``~`` included)."""
    try:
        return Path(str(given)).expanduser().is_absolute()
    except (TypeError, ValueError, RuntimeError):
        return False


def _same_file(one: Path, other: Path) -> bool:
    """Whether two paths name one file: the same text, or the same file on
    disk (a name in other letters where the disk ignores their case)."""
    if one == other:
        return True
    import os

    try:
        return os.path.samefile(one, other)
    except (OSError, ValueError):
        return False


#: Said where no home folder can be found: nothing is sent from any folder.
_NO_HOME = ("No home folder was found, so nothing is sent from any folder. Set `HOME` "
            "to your home folder and start the GUI again.")


def _through_keys(name: Any, root: Path | None = None) -> bool:
    """Whether a path as named goes through a place keys and credentials are
    kept (``.ssh`` and the like), by a part of its name or by where a link
    on the way leads (``keys -> ~/.ssh``), read from ``root`` where it is
    relative."""
    import os

    from fastmdxplora.remote.inputs import _PRIVATE_PARTS, _private, own_settings_kept

    try:
        named = Path(str(name).strip()).expanduser()
        if any(part.casefold() in _PRIVATE_PARTS for part in named.parts):
            return True
        if not named.is_absolute():
            if root is None:
                return False
            named = root / named
        named = Path(os.path.normpath(named))
        # Keys only: other programs' settings are named by the refusal.
        kept = own_settings_kept()
        if _private(named, kept):
            return True
        # Each folder on the way, where its links lead.
        for depth in range(2, len(named.parts) + 1):
            if _private(Path(*named.parts[:depth]).resolve(), kept):
                return True
    except (OSError, RuntimeError, ValueError):
        return False
    return False


#: A config read to tell a study's from a file of settings: at most this.
_CONFIG_READ_BYTES = 64 * 1024


def _a_study_config(file: Path) -> bool:
    """Whether ``file`` holds a study's config: a mapping that names its
    `system` or `systems`, read only up to :data:`_CONFIG_READ_BYTES`."""
    import os
    import stat

    import yaml

    try:
        # Opened without waiting, and judged as opened: a pipe put in its
        # place is never read, and no more than the most is read.
        handle = os.open(file, os.O_RDONLY | getattr(os, "O_NONBLOCK", 0))
        with os.fdopen(handle, "rb") as opened:
            found = os.fstat(opened.fileno())
            if not stat.S_ISREG(found.st_mode) or found.st_size > _CONFIG_READ_BYTES:
                return False
            text = opened.read(_CONFIG_READ_BYTES + 1)
        if len(text) > _CONFIG_READ_BYTES:
            return False
        raw = yaml.load(text.decode("utf-8"),
                        Loader=getattr(yaml, "CSafeLoader", yaml.SafeLoader))
    except (OSError, ValueError, yaml.YAMLError):
        return False
    from fastmdxplora.remote.inputs import names_a_study

    return names_a_study(raw)


def _as_typed(path: Path) -> str | None:
    """``path`` quoted for a POSIX shell, as a command the page shows, or
    None where no quoting is safe to show: on Windows, whose shells quote
    otherwise, and for a name holding a backslash, a quote (curly ones are
    PowerShell's quotes too) or characters that are not printable."""
    import os
    import shlex

    text = str(path)
    # A backslash or a single quote is read otherwise by some shells (fish),
    # and `!` inside single quotes by csh.
    if os.name == "nt" or not text.isprintable() or any(
            mark in text for mark in "\\'!\u2018\u2019\u201a\u201b\u201c\u201d"):
        return None
    return shlex.quote(text)


def _a_folder(path: Path) -> bool:
    """Whether ``path`` is a folder, False where it cannot be looked at."""
    try:
        return path.is_dir()
    except (OSError, ValueError):
        return False


def _said(exc: BaseException) -> dict[str, Any]:
    found = refusal_of(exc)
    return _refused(found.message, found.code)


class RemoteDesk:
    """The remote routes of one GUI server, held to its workspace: the
    folder the window was opened on and, given the window's runtime, the
    folder it puts new studies in, unless that is the home folder or one
    above it, which nothing is sent from."""

    def __init__(self, where: Any, *, hosted: bool = False) -> None:
        self.runtime = where if hasattr(where, "workspace_root") else None
        given = ([where.workspace_root, where.exploration_root] if self.runtime is not None
                 else [where])
        from fastmdxplora.remote.inputs import home_or_above

        found = [Path(root).expanduser().resolve() for root in given]
        # Left out as `gather_inputs` reads it (a disk's top, a name the
        # disk reads as home's, home not found), not by the path's text.
        self.roots = list(dict.fromkeys(
            found[:1] + [root for root in found[1:]
                         if not home_or_above(root, home_unknown=True)]))
        self.root = self.roots[0]
        self.hosted = hosted
        self._plans: dict[str, _Plan] = {}
        self._lock = threading.Lock()
        # The configs written here from the Config Builder's form, each with
        # the sha256 it was written as; those a send has begun with, kept
        # for good; and how many plans or sends are reading each now.
        self._made: dict[Path, str] = {}
        self._sent: set[Path] = set()
        self._reading: dict[Path, int] = {}

    # -- paths --------------------------------------------------------------
    def inside(self, given: Any) -> Path | None:
        """A path as named, read inside the workspace, or None: a relative
        path is the first folder's; an absolute one may be in either."""
        from fastmdxplora.mcp.workspace import Workspace

        try:
            absolute = Path(str(given)).expanduser().is_absolute()
        except (TypeError, ValueError, RuntimeError):
            return None
        for root in self.roots if absolute else self.roots[:1]:
            found = Workspace(root).inside(given)
            if found is not None:
                return found
        return None

    def _new_studies_folder(self) -> Path:
        """Where the Config Builder's studies go: the window's folder for
        them, else the folder opened."""
        return (Path(self.runtime.exploration_root).expanduser().resolve()
                if self.runtime is not None else self.root)

    def _not_built_here(self) -> str | None:
        """Why a study built in the form cannot be sent from this window,
        or None: its config would be saved in the home folder or one above
        it, which nothing is sent from. Said with the folder to open."""
        from fastmdxplora.remote.inputs import home_folders, home_or_above

        base = self._new_studies_folder()
        if not home_or_above(base, home_unknown=True):
            return None
        homes = home_folders()
        if not homes:
            # Nothing is sent where home cannot be told, a config file
            # included, so the note offers none.
            return (f"No home folder was found, so FastMDXplora cannot tell whether {base}, "
                    "where studies built in the form would be saved, is your home folder. "
                    "Nothing is sent from the GUI. Set `HOME` to your home folder and "
                    "start the GUI again.")
        still = " A config file can still be sent with A config I have."
        is_home = any(_same_file(base, home) for home in homes)
        where = ("your home folder" if is_home
                 else "the top of the disk" if base == Path(base.anchor)
                 else "a folder above your home folder")
        opened = self.root
        # The GUI puts new studies in the folder above the one it opened.
        if not home_or_above(opened, home_unknown=True):
            typed = _as_typed(opened / "first")
            instead = ("To send one, start the GUI one folder deeper, for example "
                       + (f"`fastmdx gui --output {typed}`." if typed else
                          f"with `fastmdx gui --output` naming a new folder in {opened}."))
        else:
            # The home it was opened on or in, else the account's own.
            home = next((home for home in homes if opened == home), None) or next(
                (home for home in homes if _a_folder(home)), homes[0])
            typed = _as_typed(home / "md" / "first")
            instead = ("To send one, start the GUI on a folder of its own inside your home "
                       "folder, for example "
                       + (f"`fastmdx gui --output {typed}`." if typed else
                          f"with `fastmdx gui --output` naming a new folder in {home}."))
        return (f"Studies built in the form would be saved in {base}, {where}, and "
                f"nothing is sent from there. {instead}{still}")

    def _outside(self, given: Any) -> tuple[Path, str] | None:
        """The file ``given`` names, read as a plan reads it (a relative name
        in the folder opened), and how it can be sent, where it is a config
        the routes do not reach that a refusal may say how to send
        (:meth:`_how_to_send`); else None."""
        import os

        try:
            named = Path(str(given).strip()).expanduser()
            named = Path(os.path.normpath(named if named.is_absolute() else self.root / named))
            found = named.resolve()
            if self.inside(str(found)) is None:
                how = self._how_to_send(found, named=named)
                if how is not None:
                    # Said by its name, not by where a link in it leads.
                    return named, how
        except StudyError:
            raise  # where settings are kept cannot be told: said as it is
        except (OSError, RuntimeError, ValueError):
            pass
        return None

    @staticmethod
    def _kept_apart(file: Path) -> bool:
        """Whether ``file`` is, or is in, or beside, a place keys or settings
        may be kept, and so is told no way to be sent, whatever it holds: a
        name starting with a dot; a place :func:`kept_place_of` finds (by
        name or where links lead); a folder whose name starts with a dot or
        is ``AppData`` anywhere above it; or a folder holding a place of
        keys, ``.config`` or ``.local``, or, for a home folder, ``Library``."""
        from fastmdxplora.remote.inputs import _PRIVATE_PARTS, home_folders, kept_place_of

        if file.name.startswith(".") or file.name.casefold() in _PRIVATE_PARTS:
            return True
        if kept_place_of(file, a_file=True) is not None:
            return True
        folder = file.parent
        for place in (folder, folder.resolve()):
            if any(part.startswith(".") or part.casefold() == "appdata"
                   for part in place.parts[1:]):
                return True
        # Its folder would be opened, and its pages list what is in it: each
        # name looked for, not the folder listed, which may be large.
        beside = [*_PRIVATE_PARTS, ".config", ".local"]
        if any(_same_file(folder, home) for home in home_folders()):
            beside.append("Library")
        return any((folder / name).exists() for name in beside)

    def _how_to_send(self, file: Path, started: Path | None = None, *,
                     named: Path | None = None) -> str | None:
        """How a config the routes do not reach can be sent, by ways that
        work if followed, or None where no way is said: a file that is not a
        ``.yml`` config, or one kept apart (:meth:`_kept_apart`), where it
        is and as it was named (a link in a folder of keys). The ways: a GUI
        opened on its folder, where that folder may travel, or a copy in a
        new folder inside the folder opened, where that may."""
        from fastmdxplora.remote.inputs import home_folders, home_or_above

        if not home_folders():
            # Said of every name alike, there or not: nothing is sent.
            return _NO_HOME
        try:
            if not file.is_file() or file.suffix.lower() not in (".yml", ".yaml"):
                return None
            # A study's config, by what it holds, before anything else: a
            # file of keys or settings (a kubeconfig, a token list) is none,
            # wherever it is kept and however links lead to it.
            if not _a_study_config(file):
                return None
            if any(self._kept_apart(place) for place in {file, named or file}):
                return None
            # Its folder as named, not where a link in the name leads.
            folder = (named or file).parent
            ways = []
            if not home_or_above(file.parent, home_unknown=True):
                if started is not None and _same_file(started, file.parent):
                    ways.append("start the GUI again in the config's folder without "
                                "`--output`, so it opens there")
                else:
                    typed = _as_typed(folder)
                    ways.append(f"open the GUI on the config's folder (`cd {typed}`, then "
                                "`fastmdx gui` with no `--output`)" if typed else
                                f"open the GUI on the config's folder, {folder} "
                                "(`fastmdx gui` started there, with no `--output`)")
            if self._copy_way():
                # A new folder, so nothing already there is written over.
                ways.append(f"copy the config, with the files it names, into a new folder "
                            f"inside {self.root} and give its full path there")
        except StudyError:
            raise  # where settings are kept cannot be told: said as it is
        except (OSError, RuntimeError, ValueError):
            return None
        if not ways:
            return None
        said = ", or ".join(ways) + "."
        return said[0].upper() + said[1:]

    def _copy_way(self) -> str:
        """The way any config kept elsewhere can be sent, said alike of
        every name a refusal gives no way of its own, or nothing where a new
        folder in the folder opened would be refused in turn: above the home
        folder or at a disk's top (one in the home folder itself is not)."""
        from fastmdxplora.remote.inputs import home_folders, home_or_above

        if home_or_above(self.root, home_unknown=True) and not any(
                _same_file(self.root, home) for home in home_folders()):
            return ""
        return (" To send a config kept elsewhere, copy it, with the files it names, "
                f"into a new folder inside {self.root} and give its full path there.")

    def _workspace(self) -> str:
        """The folders the routes reach, as a refusal names them."""
        if len(self.roots) == 1:
            return f"the workspace ({self.root})"
        return f"the workspace ({self.root}, and {self.roots[1]} where new studies go)"

    def shown(self, path: Any) -> str:
        """A path as the page says it: relative to the folder relative names
        are read from, else in full, so what is shown names what it says if
        it is given back."""
        from fastmdxplora.mcp.workspace import Workspace

        try:
            Path(str(path)).resolve().relative_to(self.root)
        except (ValueError, OSError, RuntimeError):
            return str(path)
        return Workspace(self.root).shown(path)

    def _results_name(self, output: Path) -> str:
        """A results folder as the form's Results folder would name it: a bare
        name where Run on this machine reads one as that folder, else in
        full."""
        if self.runtime is not None and output.parent == self._new_studies_folder():
            from fastmdxplora.gui.exploration import _slug

            if _slug(output.name) == output.name:
                return output.name
        return self.shown(output)

    def _jobs(self) -> list[Any]:
        """The jobs whose results come back into the workspace, each once."""
        from fastmdxplora.remote import api

        found: dict[str, Any] = {}
        for root in self.roots:
            for job in api.jobs(under=root):
                found.setdefault(job.name, job)
        return list(found.values())

    # -- the routes -----------------------------------------------------------
    def get(self, path: str, query: dict[str, list[str]]) -> dict[str, Any] | None:
        """The answer to a GET of ``path``, or None where it is not one of
        these routes."""
        if path not in _GETS:
            return None
        if self.hosted:
            return dict(_HOSTED)
        name = (query.get("job") or [""])[0]
        if path == "/api/remote/machines":
            return self.machines()
        if path == "/api/remote/job":
            return self.job(name)
        if path == "/api/remote/fetch-sizes":
            return self.fetch_sizes(name)
        return None

    def post(self, path: str, body: dict[str, Any]) -> dict[str, Any] | None:
        """The answer to a POST of ``path``, or None where it is not one of
        these routes."""
        if path not in _POSTS:
            return None
        if self.hosted:
            return dict(_HOSTED)
        if path == "/api/remote/plan":
            return self.plan(body.get("config"), body.get("machine"),
                             state=body.get("state"))
        if path == "/api/remote/send":
            return self.send(body.get("plan"))
        if path == "/api/remote/forget":
            return self.forget(body.get("plan"))
        if path == "/api/remote/fetch":
            return self.fetch(body.get("job"), bool(body.get("with_trajectory")),
                              body.get("bringing"))
        if path == "/api/remote/cancel":
            return self.cancel(body.get("job"))
        return None

    # -- machines and jobs ------------------------------------------------------
    def machines(self) -> dict[str, Any]:
        """The machines as recorded at a terminal, and the jobs whose results
        come back into this workspace. Nothing is asked of any machine."""
        from fastmdxplora.remote import api

        try:
            targets = [{"name": t.name, "kind": t.kind, "ready": t.ready,
                        "summary": t.summary, "gpus": list(t.gpus),
                        "inspected_at": t.inspected_at}
                       for t in api.machines()]
            jobs = [self._job_view(job) for job in self._jobs()]
        except Exception as exc:  # noqa: BLE001 - a refusal, said as one
            return _said(exc)
        return {"ok": True, "machines": targets, "jobs": jobs,
                "not_built_here": self._not_built_here()}

    def _job_here(self, given: Any):
        name = str(given or "")
        for job in self._jobs():
            if job.name == name:
                return job
        return None

    def _job_view(self, job: Any) -> dict[str, Any]:
        return {"name": job.name, "machine": job.machine, "state": job.state,
                "detail": job.detail, "scheduler": job.scheduler,
                "submitted_at": job.submitted_at, "fetched_at": job.fetched_at,
                "results": self.shown(job.local_output), "remote_dir": job.remote_dir,
                "log_tail": list(job.extra.get("log_tail") or [])[-8:]}

    def job(self, name: Any) -> dict[str, Any]:
        """A job's state, asked of its machine at most every 30 s."""
        from fastmdxplora.remote import api

        job = self._job_here(name)
        if job is None:
            return _refused(f"No job called {name!r} sends its results into this "
                            "workspace.", "remote.job.unknown")
        try:
            job = api.status(job.name, max_age_s=api.STATUS_KEPT_S)
        except Exception as exc:  # noqa: BLE001 - a refusal, said as one
            return _said(exc)
        return {"ok": True, "job": self._job_view(job)}

    # -- sending ----------------------------------------------------------------
    def _planned(self, config: Path, machine: str, output: Path):
        """The send planned, refused (as a page answer) where it may not go:
        an input outside the workspace, no room, a name or folder taken."""
        from fastmdxplora.remote import api
        from fastmdxplora.remote.jobs import job_names

        sending = api.plan_send(config, machine, output=output)
        outside = [source for source in sending.inputs.files.values()
                   if self.inside(source) is None]
        if outside:
            return None, _refused(
                f"{outside[0]} would be sent with the study, and it is outside "
                f"{self._workspace()}; nothing outside it is sent from here. Copy it "
                "into the config's folder and name it there.", "remote.input.outside")
        if sending.no_room:
            return None, _refused(sending.no_room, "remote.machine.no_room")
        if sending.job_name in job_names():
            return None, _refused(
                f"A job called {sending.job_name} was sent from this computer before, and "
                "a job's name is its results folder's. Set `output` in the config to a "
                "new folder name and plan it again.", "environment.path.exists")
        if output in self.roots or (output.exists() and (
                not output.is_dir() or any(output.iterdir()))):
            return None, _refused(
                f"{self.shown(output)} is in use already, and results are never written "
                "over anything. Set `output` in the config to a new folder name and plan "
                "it again.", "environment.path.exists")
        return sending, None

    def _what(self, config: Any) -> tuple[Path | None, Path | None, dict[str, Any] | None]:
        """The config file and its results folder, both in the workspace."""
        import yaml

        # A name through a place of keys reads as any name that is no
        # config, wherever its links lead.
        through = bool(config) and _through_keys(config, self.root)
        file = self.inside(config) if config and not through else None
        if file is None and config and not through:
            outside = self._outside(config)
            if outside is not None:
                found, how = outside
                named = (f"{found} is" if str(found) == str(config).strip()
                         else f"{config} is {found},")
                return None, None, _refused(
                    f"{named} outside {self._workspace()}, and nothing outside it is "
                    f"sent from here. {how}", "remote.input.outside")
        if file is None or not file.is_file() or file.suffix.lower() not in (".yml", ".yaml"):
            # Said alike of every such name, there or not, inside or not.
            return None, None, _refused(
                f"{config} is not a config file in " + (
                    self._workspace() if _absolute(config) else str(self.root)) + "."
                + self._copy_way())
        import os

        from fastmdxplora.remote.inputs import kept_place_of, kept_place_said

        def as_given(name: Any) -> Path:
            named = Path(str(name).strip()).expanduser()
            return Path(os.path.normpath(named if named.is_absolute() else self.root / named))

        # In a window on the home folder its places of settings are inside:
        # judged as named and where links lead, before the file is read.
        place = kept_place_of(as_given(config), a_file=True) or kept_place_of(
            file, a_file=True)
        if place is not None:
            # Inside the workspace, where the page may look: a file that is
            # no study's config is called a file.
            return None, None, _refused(
                kept_place_said(str(config).strip(), place, a_study=_a_study_config(file)),
                "remote.input.outside")
        try:
            raw = yaml.safe_load(file.read_text(encoding="utf-8"))
        except (OSError, ValueError, yaml.YAMLError) as exc:
            # What was wrong, without the file's own lines.
            why = getattr(exc, "problem", None) or type(exc).__name__
            return None, None, _refused(f"{self.shown(file)} could not be read: {why}. "
                                        "Correct the file and plan it again.")
        requested = str((raw or {}).get("output") or file.with_suffix("")) if isinstance(
            raw, dict) else str(file.with_suffix(""))
        output = self.inside(requested)
        if output is None:
            return None, None, _refused(f"The results folder {requested} is outside "
                                        f"{self._workspace()}. Set `output` in the config "
                                        f"to a new folder inside it.", "remote.input.outside")
        place = kept_place_of(as_given(requested)) or kept_place_of(output)
        if place is not None:
            return None, None, _refused(
                kept_place_said(requested, place, results=True, in_gui=True),
                "remote.input.outside")
        return file, output, None

    def _saved(self, state: dict[str, Any]
               ) -> tuple[Path | None, dict[str, Any] | None, bool]:
        """The Config Builder's study saved as a config file, named after the
        results folder it names (``<results>.yml``, or ``-2``, ``-3`` and so
        on beside an earlier one), in the workspace's folder: the folder new
        studies go in, given the window's runtime, else the folder opened.
        A file the study names by a relative path is refused where **Run on
        this machine** would read another, or none (it starts the run in the
        folder new studies go in, the one above the folder opened as
        ``fastmdx gui`` sets it), so what travels is never a different file
        under the same name; its full path names it either way. A config
        file the person has is sent as ``fastmdx remote send -c`` sends it,
        with the files beside it, each shown in the plan. The same
        text planned again is the same file; a link, or a file that is not
        the same text, is passed over and never read or written. Also
        whether the file was written now (and so is taken back where the
        plan is refused)."""
        from fastmdxplora.gui.config_builder import config_yaml
        from fastmdxplora.naming import default_output_name, system_of

        base = self._new_studies_folder()
        why = self._not_built_here()
        if why is not None:
            # Nothing in the home folder itself travels, and a config is not
            # left there to say so.
            return None, _refused(why), False
        if base not in self.roots:
            # Only into a folder the routes reach, however the window's
            # folder for new studies came to be another since.
            return None, _refused(f"Studies built in the form would be saved in {base}, which "
                                  "this window does not send from. Start the GUI again "
                                  "to send one."), False
        source = dict(state)
        requested = str(source.get("output") or "").strip()
        if not requested:
            requested = default_output_name(system_of(source))
        try:
            output = (self.runtime._output_folder(requested) if self.runtime is not None
                      else self.inside(requested))
        except (OSError, ValueError, RuntimeError):
            output = None
        if output is None or self.inside(str(output)) is None or output in self.roots:
            return None, _refused(f"The results folder {requested} must be a new folder "
                                  f"inside {self._workspace()}."), False
        source["output"] = str(output)
        try:
            built = config_yaml(source, full=bool(source.get("full")))
        except Exception as exc:  # noqa: BLE001 - a refusal, said as one
            return None, _said(exc), False
        if not built.get("ok"):
            return None, _refused(str(built.get("error") or "The study could not be "
                                      "written as a config."), str(built.get("code") or "")), False
        text = str(built["yaml"])
        run_reads = (Path(self.runtime.exploration_root).expanduser().resolve()
                     if self.runtime is not None else self.root.parent)
        if run_reads != base:
            named = _read_otherwise(text, base, run_reads)
            if named:
                from fastmdxplora.remote.inputs import _existing

                none_here = ", where there is none," if _existing(
                    named[0], run_reads) is None else ""
                return None, _refused(
                    f"Run on this machine would read {named[0]!r} from {run_reads}"
                    f"{none_here} but a send reads it from {base}. Give its full path "
                    "in the form (Browse beside the field gives one), so the file that "
                    "travels is the one a run here would read.",
                    "remote.input.outside"), False
        for n in range(1, 1000):
            target = base / (f"{output.name}.yml" if n == 1 else f"{output.name}-{n}.yml")
            try:
                if target.is_symlink() or (target.exists() and (
                        self.inside(str(target)) is None or not target.is_file())):
                    continue  # never read, nor written through
                if target.is_file():
                    try:
                        if target.read_text(encoding="utf-8") == text:
                            return target, None, False
                    except ValueError:  # not text: another file, kept as it is
                        pass
                    continue
                with target.open("x", encoding="utf-8") as out:
                    out.write(text)
            except FileExistsError:
                continue
            except (OSError, ValueError) as exc:
                return None, _refused(f"The config could not be saved in "
                                      f"{self.shown(base)}: {exc}"), False
            return target, None, True
        return None, _refused("Too many configs of that name; give another results "
                              "folder."), False

    def plan(self, config: Any, machine: Any, *, state: Any = None) -> dict[str, Any]:
        """What a send of ``config`` to ``machine`` would do, as the page
        shows it before the person agrees; with a token the send is made by.
        Given the Config Builder's ``state`` in place of a file, the study is
        saved as one first (:meth:`_saved`). A refusal on the way (where
        settings are kept cannot be told, say) is said as one."""
        try:
            return self._plan_of(config, machine, state=state)
        except StudyError as exc:
            return _said(exc)

    def _plan_of(self, config: Any, machine: Any, *, state: Any = None) -> dict[str, Any]:
        from fastmdxplora.remote.inputs import home_folders

        if isinstance(state, dict) and state.get("review_after_setup") is True:
            return _refused("Review after setup runs locally on this machine; it cannot be sent elsewhere.",
                            "config.option.inapplicable")
        if not home_folders():
            # Nothing is sent from any folder, so nothing more is looked at
            # (a machine's record is found through the home folder too).
            return _refused(_NO_HOME, "remote.input.outside")
        if not isinstance(state, dict) and not str(config or "").strip():
            return _refused("No config file was named. Give the full path of the config "
                            "to send.")
        if not config and isinstance(state, dict):
            from fastmdxplora.remote.machines import load_machine

            # A machine never inspected is said before anything is saved.
            try:
                load_machine(str(machine or ""))
            except Exception as exc:  # noqa: BLE001 - a refusal, said as one
                return _said(exc)
            with self._lock:
                # Saved and counted as read at once, so a plan let go in
                # another tab never takes back the file this one found.
                saved, refused, made = self._saved(state)
                if refused is None:
                    if made:
                        self._made[saved] = _sha_of(saved)
                    self._reading[saved] = self._reading.get(saved, 0) + 1
            if refused is not None:
                return refused
            try:
                said = self._plan(str(saved), machine)
            finally:
                self._done_reading(saved)
            if not said.get("ok"):
                # Refused: its config is taken back, so planned again it is
                # written again, and a blank results folder leaves none.
                self._release(saved)
            return said
        refused = self._read_elsewhere_here(config)
        if refused is not None:
            return refused
        # A config the person has may be one the form wrote: held while it
        # is planned, as the form's own plan holds it.
        named = self.inside(config) if config else None
        if named is None:
            return self._plan(config, machine)
        with self._lock:
            self._reading[named] = self._reading.get(named, 0) + 1
        try:
            said = self._plan(config, machine)
        finally:
            self._done_reading(named)
        if not said.get("ok"):
            # A form's config let go while this plan held it is taken back
            # now, as it would have been then.
            self._release(named)
        return said

    def _read_elsewhere_here(self, config: Any) -> dict[str, Any] | None:
        """A config file named by a relative path is refused where the
        window checks and runs another, or none: **A config I have** reads
        such a name from the folder the GUI was started in, and a plan of
        it from the folder it was opened on, so what is planned is always
        the file checked. Refused too where the folder the GUI was started
        in is gone, since the file checked can no longer be told."""
        if self.runtime is None or not config:
            return None
        text = str(config).strip()
        try:
            given = Path(text)
            if text.startswith("~") or given.is_absolute():
                return None
        except (RuntimeError, ValueError):
            return None
        code = "remote.config.read_elsewhere"
        full = " Give its full path, so the file planned is the one checked."
        try:
            started = Path.cwd().resolve()
        except (OSError, RuntimeError):
            return _refused(f"{text} is checked from the folder the GUI was started in, "
                            "and that folder is gone, so which file was checked is not "
                            "known." + full, code, given=text)
        try:
            checked = (started / given).resolve()
        except ValueError:
            return None  # no file's name: refused by the plan as no config file
        except (OSError, RuntimeError) as exc:
            # A loop of links on the way, say: the file checked is not known.
            return _refused(f"{text} cannot be followed from the folder the GUI was "
                            f"started in ({started}): {exc}. So which file was checked "
                            "is not known." + full, code, given=text)
        import os

        # As typed, read from each folder by the name's text alone: where a
        # link outside the routes leads is not said.
        typed = Path(os.path.normpath(started / given))
        there = Path(os.path.normpath(self.root / given))
        planned = self.inside(text)
        # The same file in the same folder: the files it names by relative
        # paths are then the same files too (a hard link elsewhere is not).
        if planned is not None and _same_file(checked, planned) and _same_file(
                checked.parent, planned.parent):
            return None

        # A name through a place of keys is not followed, even into the
        # routes: what is in such a place is never told.
        inside = (self.inside(str(checked)) is not None
                  and not _through_keys(text, started))  # what is in it is not told
        try:
            # A link crossed before a `..` (`.ssh/lnk/../x.yml`), or a name
            # through a folder kept apart by its dot: no way is offered.
            crossed = checked != typed.resolve()
        except (OSError, RuntimeError, ValueError):
            crossed = True
        dotted = any(part.startswith(".") and part not in (".", "..") for part in given.parts)
        through = _through_keys(text, started)
        how = (None if inside or crossed or dotted or through
               else self._how_to_send(checked, started, named=typed))
        # Where a link leads is said only inside the routes.
        shown = checked if inside else typed
        read = (f"Checked from the folder the GUI was started in ({started}), {text} is "
                f"{shown}; ")
        try:
            # Inside the routes, where the page may look anyway.
            there_now = planned is not None and planned.exists()
        except (OSError, ValueError):
            there_now = False
        if planned is None:
            # A plan reads a relative name only inside the folder opened.
            said = (f"a plan reads a relative name only inside the folder the GUI was "
                    f"opened on ({self.root}), and {text} leads out of it.")
            said = (said[0].upper() + said[1:]) if there == typed else read + said
        elif there_now:
            said = (read + f"a plan reads it from the folder the GUI was opened on "
                    f"({self.root}), as {planned}.")
        else:
            said = (read + f"a plan reads it from the folder the GUI was opened on "
                    f"({self.root}), where there is no such file.")
        if inside:
            try:
                reached = checked.is_file()
            except (OSError, ValueError):
                reached = False
            if reached:
                return _refused(said + full, code, given=text)
        else:
            if how is not None:
                # A config its full path would not reach either: how to send it.
                it = ("The file checked" if said.startswith("Checked")
                      else f"The file {text} names")
                return _refused(f"{said} {it} is outside {self._workspace()}, "
                                "so its full path would be refused too: nothing outside "
                                f"the workspace is sent from here. {how}", code, given=text)
        # Said alike of any other file and of none, so nothing more is told.
        copy = "" if inside else self._copy_way()
        if copy:
            # The full path of a file outside would be refused too.
            return _refused(said + " Give the full path of a config in the workspace, or,"
                            + copy.replace(" To send", " to send", 1), code, given=text)
        return _refused(said + " Give the full path of the config to send.", code,
                        given=text)

    def _plan(self, config: Any, machine: Any) -> dict[str, Any]:
        from fastmdxplora.remote.send import room_said, sent_digest, travelling

        file, output, refused = self._what(config)
        if refused is not None:
            return refused
        name = str(machine or "")
        try:
            sending, refused = self._planned(file, name, output)
        except Exception as exc:  # noqa: BLE001 - a refusal, said as one
            return _said(exc)
        if refused is not None:
            return refused
        prints = travelling(sending)
        token = secrets.token_hex(16)
        with self._lock:
            dropped = self._forget_old()
            while len(self._plans) >= _MOST_PLANS:
                dropped.append(self._plans.pop(
                    min(self._plans, key=lambda t: self._plans[t].made_at)).config)
            self._plans[token] = _Plan(file, name, output, sent_digest(sending, prints),
                                        time.monotonic())
        for gone in dropped:
            self._release(gone)
        return {
            "ok": True, "plan": token, "kept_s": PLAN_KEPT_S,
            "config": self.shown(file), "machine": name, "job": sending.job_name,
            "runs_in": sending.installation.path,
            "scheduler": "SLURM" if sending.scheduler == "slurm" else "a detached process",
            "folder": sending.remote_dir, "results": self.shown(output),
            "results_path": str(output), "results_name": self._results_name(output),
            "travels": [{"name": f"inputs/{travelled}", "from": self.shown(source),
                         "bytes": prints[travelled][0]}
                        for travelled, source in sending.inputs.files.items()],
            "fetched_there": list(sending.inputs.fetched),
            "room": room_said(sending), "notes": list(sending.notes),
            "script": sending.script,
        }

    def forget(self, token: Any) -> dict[str, Any]:
        """A plan the page let go (Not now, or the study changed): never
        sent, and the config written for it taken back (:meth:`_release`)."""
        with self._lock:
            plan = self._plans.pop(str(token or ""), None)
        if plan is not None:
            self._release(plan.config)
        return {"ok": True}

    def _done_reading(self, config: Path) -> None:
        with self._lock:
            left = self._reading.get(config, 0) - 1
            if left > 0:
                self._reading[config] = left
            else:
                self._reading.pop(config, None)

    def _release(self, config: Path) -> None:
        """A config written here from the form, taken back once nothing
        holds it: no plan of it kept, none being made or sent, no send ever
        begun with it, and the file as it was written."""
        with self._lock:
            if (config not in self._made or config in self._sent
                    or self._reading.get(config, 0)
                    or any(p.config == config for p in self._plans.values())):
                return
            _taken_back(config, self._made.pop(config))

    def _forget_old(self) -> list[Path]:
        """Plans kept past their ten minutes dropped (under the lock); their
        configs, for :meth:`_release` once the lock is let go."""
        now = time.monotonic()
        return [self._plans.pop(token).config for token in
                [t for t, p in self._plans.items() if now - p.made_at > PLAN_KEPT_S]]

    def send(self, token: Any) -> dict[str, Any]:
        """Send the plan the person agreed to, once, where what would travel
        now is what they were shown."""
        from fastmdxplora.remote import api
        from fastmdxplora.remote.send import sent_digest

        with self._lock:
            dropped = self._forget_old()
            plan = self._plans.pop(str(token or ""), None)
            if plan is not None:
                self._reading[plan.config] = self._reading.get(plan.config, 0) + 1
        for gone in dropped:
            self._release(gone)
        if plan is None:
            return _refused("That plan is not one this page was shown in the last ten "
                            "minutes, or it was sent already. Plan the send again.",
                            "remote.send.unconfirmed")
        try:
            sending, refused = self._planned(plan.config, plan.machine, plan.output)
        except Exception as exc:  # noqa: BLE001 - a refusal, said as one
            sending, refused = None, _said(exc)
        if refused is None and sent_digest(sending) != plan.digest:
            refused = _refused("What would be sent changed after the plan was shown (the "
                               "config, a file that travels, or the machine). Plan the "
                               "send again and look at it.", "remote.send.unconfirmed")
        with self._lock:
            if refused is None:
                # From here the config is a sent job's, or may be: kept.
                self._sent.add(plan.config)
        self._done_reading(plan.config)
        if refused is not None:
            # Nothing went: the plan is let go as Not now lets it go.
            self._release(plan.config)
            return {**refused, "sent_nothing": True}
        try:
            job = api.send_planned(sending)
        except Exception as exc:  # noqa: BLE001 - a refusal, said as one
            return _said(exc)
        return {"ok": True, "job": self._job_view(job)}

    # -- fetching and stopping ------------------------------------------------
    def fetch_sizes(self, name: Any) -> dict[str, Any]:
        """What a fetch of an ended job would bring, each part's size."""
        from fastmdxplora.remote import api

        job = self._job_here(name)
        if job is None:
            return _refused(f"No job called {name!r} sends its results into this "
                            "workspace.", "remote.job.unknown")
        try:
            job = api.status(job.name, max_age_s=api.STATUS_KEPT_S)
            if job.state in ("ready", "running"):
                return _refused(f"{job.name} is still {job.state}; it is fetched once it "
                                "has ended.", "remote.job.unfinished")
            sizes = api.fetch_sizes(job.name, max_age_s=api.STATUS_KEPT_S)
        except Exception as exc:  # noqa: BLE001 - a refusal, said as one
            return _said(exc)
        return {"ok": True, "job": self._job_view(job), "results_bytes": sizes.results,
                "trajectory_bytes": sizes.trajectory,
                "trajectory_files": sizes.trajectory_files,
                "run_written": sizes.run_written,
                "bringing": {"without": sizes.bringing(False), "with": sizes.bringing(True)}}

    def fetch(self, name: Any, with_trajectory: bool, bringing: Any) -> dict[str, Any]:
        """Fetch a job's results, where what it would bring, asked of the
        machine now, is still what the person was shown (``bringing``, in
        bytes). The copy then brings no file larger than that, and a margin."""
        from fastmdxplora.remote import api

        job = self._job_here(name)
        if job is None:
            return _refused(f"No job called {name!r} sends its results into this "
                            "workspace.", "remote.job.unknown")
        try:
            sizes = api.fetch_sizes(job.name)
        except Exception as exc:  # noqa: BLE001 - a refusal, said as one
            return _said(exc)
        now = sizes.bringing(with_trajectory)
        if not (isinstance(bringing, int) and not isinstance(bringing, bool)
                and bringing == now):
            return _refused("The fetch would bring another size than the one shown. Look "
                            "at the sizes again.", "remote.fetch.unconfirmed",
                            bringing=now)
        try:
            # No file larger than the whole was said to be, and a margin.
            job, warnings = api.fetch(job.name, with_trajectory=with_trajectory,
                                      most_bytes=now + now // 10 + 1_000_000)
        except Exception as exc:  # noqa: BLE001 - a refusal, said as one
            return _said(exc)
        return {"ok": True, "job": self._job_view(job), "warnings": list(warnings)}

    def cancel(self, name: Any) -> dict[str, Any]:
        """Stop a job, as `fastmdx remote cancel` does (the page asks first)."""
        from fastmdxplora.remote import api
        from fastmdxplora.remote.send import left_stop_done

        job = self._job_here(name)
        if job is None:
            return _refused(f"No job called {name!r} sends its results into this "
                            "workspace.", "remote.job.unknown")
        asked_at = time.time()
        try:
            job = api.cancel(job.name)
        except Exception as exc:  # noqa: BLE001 - a refusal, said as one
            return _said(exc)
        return {"ok": True, "job": self._job_view(job),
                "stopped": job.state == "abandoned" or left_stop_done(job, asked_at) in (
                    "sent", "now")}
