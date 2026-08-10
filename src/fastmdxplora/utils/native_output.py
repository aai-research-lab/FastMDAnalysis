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
import sys
from collections.abc import Iterator


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

    saved: dict[int, int] = {}
    devnull_fd: int | None = None

    def _flush_all() -> None:
        # Flush Python streams, then the C runtime's stdio buffers. The
        # libc flush is what catches buffered C printf output on macOS.
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

    try:
        _flush_all()
        try:
            devnull_fd = os.open(
                os.devnull if into is None else os.fspath(into),
                os.O_WRONLY | (0 if into is None
                               else os.O_CREAT | os.O_APPEND))
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
