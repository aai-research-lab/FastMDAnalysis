"""Suppress raw C-level stdout/stderr from extension libraries.

Some C extensions — notably MDTraj's DCD/PDB molfile plugins — write
status messages straight to the operating-system file descriptors
(``dcdplugin) detected standard 32-bit DCD file ...``). Because the
writes happen in C, below Python, they bypass ``sys.stdout`` /
``sys.stderr`` and Python's logging entirely; only redirecting the
underlying file descriptors (1 and 2) intercepts them.

This module provides one context manager, :func:`suppress_native_output`,
used both by the trajectory loader (so real ``fastmdx`` runs stay clean)
and by the test fixtures (which build/read DCD files directly).
"""

from __future__ import annotations

import contextlib
import ctypes
import ctypes.util
import os
import re
import sys
import threading
import warnings
from collections.abc import Iterator
from typing import Any


def _load_libc() -> ctypes.CDLL | None:
    """Load libc for flushing C stdio buffers, or None if unavailable."""
    try:
        name = ctypes.util.find_library("c")
        # On some platforms find_library returns None but the default CDLL
        # (None) still resolves the C runtime; try both.
        return ctypes.CDLL(name) if name else ctypes.CDLL(None)
    except (OSError, TypeError):
        return None


_LIBC = _load_libc()


_REDIRECT_HELD = threading.Lock()
_redirect_count: dict[int, int] = {}
_redirect_saved: dict[int, int] = {}
_redirect_null: list[int] = []


def _flush_all() -> None:
    """Python's streams, then the C runtime's stdio buffers (the libc flush
    is what catches buffered C printf output on macOS)."""
    try:
        sys.stdout.flush()
        sys.stderr.flush()
    except (ValueError, OSError):
        pass
    if _LIBC is not None:
        try:
            _LIBC.fflush(None)  # fflush(NULL) flushes ALL open streams
        except OSError:
            pass


@contextlib.contextmanager
def _shared_redirect(fds: list[int]) -> Iterator[None]:
    """``fds`` sent to /dev/null for the block: the first to ask redirects,
    the last to leave puts the terminal back, the ones between only count."""
    taken: list[int] = []
    with _REDIRECT_HELD:
        _flush_all()
        if not _redirect_null:
            try:
                _redirect_null.append(os.open(os.devnull, os.O_WRONLY))
            except OSError:
                pass
        if _redirect_null:
            for fd in fds:
                if not _redirect_count.get(fd):
                    try:
                        _redirect_saved[fd] = os.dup(fd)
                        os.dup2(_redirect_null[0], fd)
                    except OSError:
                        # This descriptor isn't redirectable; skip it.
                        _redirect_saved.pop(fd, None)
                        continue
                _redirect_count[fd] = _redirect_count.get(fd, 0) + 1
                taken.append(fd)
    try:
        yield
    finally:
        with _REDIRECT_HELD:
            # Flushed while still redirected, so buffered C output goes to
            # /dev/null, then the real descriptors put back by the last out.
            _flush_all()
            for fd in taken:
                _redirect_count[fd] -= 1
                if not _redirect_count[fd]:
                    original = _redirect_saved.pop(fd)
                    try:
                        os.dup2(original, fd)
                    finally:
                        os.close(original)
            if _redirect_null and not any(_redirect_count.values()):
                os.close(_redirect_null.pop())


@contextlib.contextmanager
def suppress_native_output(
    *, stdout: bool = True, stderr: bool = True, into: "os.PathLike | str | None" = None
) -> Iterator[None]:
    """Redirect OS-level stdout/stderr away from the terminal for the block.

    With ``into``, the output is written to that file rather than discarded.
    PLUMED prints forty lines of setup at the moment its force is added --
    which atoms the collective variable is built from, the hill width, the
    pace, the bias factor, the temperature it inferred. That is provenance,
    and the wrong thing to do with it is either leave it scrolling past a
    progress bar or throw it away; it belongs in a file beside the run.

    Operates on the file descriptors themselves (``os.dup2``), so it
    catches output written by C extensions that bypass Python's stream
    objects. Crucially, it also flushes the C runtime's stdio buffers
    (``libc.fflush(NULL)``) *while the redirect is still in place*, before
    restoring the descriptors — without this, libc-buffered C output (e.g.
    MDTraj's plugin messages on macOS) leaks to the terminal after the
    descriptors are put back. Python-level ``print`` is also suppressed
    for the duration.

    Parameters
    ----------
    stdout, stderr : bool
        Which OS streams to redirect. Both default to True because the
        target chatter has been observed on either stream depending on
        platform and how the terminal wired the descriptors.

    Notes
    -----
    Best-effort: if a descriptor can't be duplicated (e.g. it's already
    closed, as can happen under some capture setups), that stream is left
    alone rather than raising.
    """
    fds: list[int] = []
    if stdout:
        fds.append(1)
    if stderr:
        fds.append(2)

    if into is None:
        # Discarded: one redirect for the whole process, however many
        # threads ask for it at once (the GUI's server reads on several).
        # Each saving and putting back its own, two that overlapped left
        # the terminal on /dev/null for good: the second saved the first's
        # redirect as the terminal.
        with _shared_redirect(fds):
            yield
        return

    saved: dict[int, int] = {}
    devnull_fd: int | None = None

    try:
        _flush_all()
        try:
            devnull_fd = os.open(os.fspath(into), os.O_WRONLY | os.O_CREAT | os.O_APPEND)
        except OSError:
            # Can't open the destination — give up quietly rather than
            # losing the run over where a log line went.
            yield
            return
        for fd in fds:
            try:
                saved[fd] = os.dup(fd)
                os.dup2(devnull_fd, fd)
            except OSError:
                # This descriptor isn't redirectable; skip it.
                saved.pop(fd, None)
        yield
    finally:
        # Flush *while still redirected* so buffered C output goes to
        # devnull, THEN restore the real descriptors.
        _flush_all()
        for fd, original in saved.items():
            try:
                os.dup2(original, fd)
            finally:
                os.close(original)
        if devnull_fd is not None:
            os.close(devnull_fd)


# MDTraj's notes on reading a structure that change no atom: a CRYST1 cell
# too small for the atoms in it (MDTraj's test is a density above 1000 atoms
# per nm^3; the 1 A placeholder NMR entries such as 1L2Y carry is one), which
# MDTraj discards, and two residues one after another under one number and
# insertion code but different names, which it keeps as two residues (a
# structure written without the insertion codes that told them apart).
_NOTES_ON_READING = (
    "Unlikely unit cell vectors detected in PDB file likely resulting from a "
    "dummy CRYST1 record",
    "WARNING: two consecutive residues with same number",
)
_NOTES_HELD = threading.Lock()


def keep_notes_on_reading_quiet() -> None:
    """Ignore MDTraj's notes on reading a structure, from now on, in this
    process.

    The two filters are added once, by their message and for MDTraj's own
    modules, and left in place: set and put back around each load
    (``warnings.catch_warnings``), they were put back by one server thread
    while another was still reading, so that thread's note was printed, and
    its own filters were left behind for good. Called again, it adds only a
    filter that is missing (a test's own ``catch_warnings`` takes them away
    when it ends). Every other warning shows as it would, and once per
    place as Python shows it.
    """
    with _NOTES_HELD:
        # A filter's message and module are patterns, or plain text where
        # one was put in the list as text (Python's own '__main__' is).
        here = {(getattr(one[1], "pattern", one[1]), getattr(one[3], "pattern", one[3]))
                for one in warnings.filters if one[0] == "ignore"}
        for note in _NOTES_ON_READING:
            if (re.escape(note), _MDTRAJ) not in here:
                warnings.filterwarnings("ignore", message=re.escape(note), module=_MDTRAJ)


_MDTRAJ = r"mdtraj\."


def topology_quietly(path: "os.PathLike | str") -> Any:
    """The topology MDTraj reads from ``path``, its notes on reading a
    structure kept out of the terminal (:func:`keep_notes_on_reading_quiet`);
    an error raises as it would."""
    import mdtraj as md

    keep_notes_on_reading_quiet()
    return md.load_topology(os.fspath(path))


def pdb_quietly(path: "os.PathLike | str") -> Any:
    """The PDB at ``path`` as MDTraj reads it, its frames with their
    topology, its notes on reading a structure kept out of the terminal. The
    cell is still discarded as MDTraj discards it."""
    import mdtraj as md

    keep_notes_on_reading_quiet()
    return md.load_pdb(os.fspath(path))
