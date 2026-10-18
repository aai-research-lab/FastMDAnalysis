"""Request-local admission for existing record readers.

Ordinary callers retain pathlib's behavior. A public record aggregate can
enter ``confined_records`` so its readers and their helpers share one
boundary without caching a tree scan or changing another request's policy.
"""

from __future__ import annotations

import errno
import fnmatch
import os
import stat as stat_module
from contextlib import contextmanager
from contextvars import ContextVar
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Iterator, TextIO

from fastmdxplora.refusals import OutsideWorkspace


@dataclass
class ReadScope:
    root: Path
    admitted: Callable[[Path], bool]
    denied: set[str] = field(default_factory=set)
    descriptor: int | None = None

    def checked(self, path: Path, *, remember: bool = True) -> Path:
        """Check both the requested address and its current destination."""
        path = path.absolute()
        relative: Path | None = None
        try:
            relative = path.relative_to(self.root)
            resolved = path.resolve()
            resolved.relative_to(self.root)
            if self.admitted(path) and self.admitted(resolved):
                return resolved
        except (OSError, RuntimeError, ValueError):
            pass
        if remember:
            # A link out of the study is named by its own place in it; only
            # an address that is not in the study at all goes unnamed.
            named = relative is not None and ".." not in relative.parts
            self.denied.add(relative.as_posix() if named else "an unavailable record")
        raise OutsideWorkspace("The record could not be confined to the allowed study files.")


_scope: ContextVar[ReadScope | None] = ContextVar("record_read_scope", default=None)


def _root_descriptor(root: Path) -> int:
    flags = os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW
    # A folder above the study need only be searchable (a shared cluster's
    # home folder is often 0711). Where the system has O_PATH, those are
    # opened without asking to read them; the study itself is opened to read.
    above = (os.O_PATH | os.O_DIRECTORY | os.O_NOFOLLOW) if hasattr(os, "O_PATH") else flags
    parts = root.parts[1:]
    descriptor = os.open(root.anchor, above if parts else flags)
    try:
        for index, part in enumerate(parts):
            final = index == len(parts) - 1
            opened = os.open(part, flags if final else above, dir_fd=descriptor)
            os.close(descriptor)
            descriptor = opened
        return descriptor
    except BaseException:
        os.close(descriptor)
        raise


@contextmanager
def confined_records(root: Path, admitted: Callable[[Path], bool]) -> Iterator[ReadScope]:
    """Apply an admission policy only to this call's record reads."""
    scope = ReadScope(Path(root).resolve(), admitted)
    # Open relative to a pinned directory on POSIX. Each resolved component
    # is opened without following links, so retargeting one after admission
    # cannot redirect an actual read outside the admitted tree.
    if os.open in os.supports_dir_fd and hasattr(os, "O_NOFOLLOW"):
        scope.descriptor = _root_descriptor(scope.root)
    token = _scope.set(scope)
    try:
        yield scope
    finally:
        _scope.reset(token)
        if scope.descriptor is not None:
            os.close(scope.descriptor)


def _open(scope: ReadScope, path: Path, *, directory: bool = False) -> int:
    resolved = scope.checked(path)
    flags = os.O_RDONLY | getattr(os, "O_NONBLOCK", 0)
    if scope.descriptor is not None:
        parts = resolved.relative_to(scope.root).parts
        descriptor = os.dup(scope.descriptor)
        try:
            for index, part in enumerate(parts):
                final = index == len(parts) - 1
                opened = os.open(part, flags | os.O_NOFOLLOW |
                                 (os.O_DIRECTORY if directory or not final else 0),
                                 dir_fd=descriptor)
                os.close(descriptor)
                descriptor = opened
            return descriptor
        except OSError as exc:
            os.close(descriptor)
            # Record a rejected link even when an inner legacy reader
            # catches the resulting OSError and supplies a fallback.
            scope.checked(path)
            if exc.errno in {errno.ELOOP, errno.ENOTDIR}:
                scope.denied.add(path.relative_to(scope.root).as_posix())
            raise
    descriptor = os.open(resolved, flags)
    try:
        # Windows does not support openat. Check the destination of the
        # opened handle itself before any bytes are read.
        import ctypes
        import msvcrt

        get_name = ctypes.windll.kernel32.GetFinalPathNameByHandleW
        get_name.argtypes = [ctypes.c_void_p, ctypes.c_wchar_p, ctypes.c_ulong,
                             ctypes.c_ulong]
        get_name.restype = ctypes.c_ulong
        handle = msvcrt.get_osfhandle(descriptor)
        length = get_name(handle, None, 0, 0)
        buffer = ctypes.create_unicode_buffer(length + 1)
        found = get_name(handle, buffer, len(buffer), 0)
        if not found or found >= len(buffer):
            scope.denied.add("an unavailable record")
            raise OutsideWorkspace("The opened record's location could not be verified.")
        name = buffer.value
        if name.startswith("\\\\?\\UNC\\"):
            name = "\\\\" + name[8:]
        elif name.startswith("\\\\?\\"):
            name = name[4:]
        scope.checked(Path(name))
        return descriptor
    except BaseException:
        os.close(descriptor)
        raise


def open_text(path: Path, *, encoding: str = "utf-8", errors: str | None = None,
              newline: str | None = None) -> TextIO:
    scope = _scope.get()
    if scope is None:
        return Path(path).open(encoding=encoding, errors=errors, newline=newline)
    descriptor = _open(scope, Path(path))
    if not stat_module.S_ISREG(os.fstat(descriptor).st_mode):
        os.close(descriptor)
        scope.denied.add("a non-regular record")
        raise OutsideWorkspace("Only regular study records may be read.")
    return os.fdopen(descriptor, encoding=encoding, errors=errors, newline=newline)


def refused() -> frozenset[str]:
    """The records this call's held read has refused so far, each named in
    the folder it is held to; none outside a held read."""
    scope = _scope.get()
    return frozenset(scope.denied) if scope is not None else frozenset()


def read_text(path: Path, *, encoding: str = "utf-8", errors: str | None = None) -> str:
    with open_text(path, encoding=encoding, errors=errors) as stream:
        return stream.read()


def can_read(path: Path) -> bool:
    """Whether a possible source is in scope, without reading or refusing it.

    This lets a reader distinguish an unknown owning-study relationship
    from a missing record. Actual reads still enforce admission afresh.
    """
    scope = _scope.get()
    if scope is None:
        return True
    try:
        scope.checked(Path(path), remember=False)
        return True
    except OutsideWorkspace:
        return False


def stat(path: Path) -> os.stat_result:
    scope = _scope.get()
    if scope is None:
        return Path(path).stat()
    if scope.descriptor is None:
        return scope.checked(Path(path)).stat()
    resolved = scope.checked(Path(path))
    if resolved == scope.root:
        return os.fstat(scope.descriptor)
    descriptor = _open(scope, resolved.parent, directory=True)
    try:
        result = os.stat(resolved.name, dir_fd=descriptor, follow_symlinks=False)
        if stat_module.S_ISLNK(result.st_mode):
            scope.denied.add(Path(path).relative_to(scope.root).as_posix())
            raise OutsideWorkspace("The record changed after its location was checked.")
        return result
    finally:
        os.close(descriptor)


def is_file(path: Path) -> bool:
    if _scope.get() is None:
        return Path(path).is_file()
    try:
        return stat_module.S_ISREG(stat(path).st_mode)
    except OSError:
        return False


def is_dir(path: Path) -> bool:
    if _scope.get() is None:
        return Path(path).is_dir()
    try:
        return stat_module.S_ISDIR(stat(path).st_mode)
    except OSError:
        return False


def exists(path: Path) -> bool:
    if _scope.get() is None:
        return Path(path).exists()
    try:
        stat(path)
        return True
    except OSError:
        return False


def iterdir(path: Path) -> Iterator[Path]:
    scope = _scope.get()
    if scope is None:
        yield from Path(path).iterdir()
        return
    path = Path(path)
    if scope.descriptor is not None:
        descriptor = _open(scope, path, directory=True)
        try:
            names = os.listdir(descriptor)
        finally:
            os.close(descriptor)
    else:
        names = os.listdir(scope.checked(path))
    for name in names:
        candidate = path / name
        # Discovery ignores private names. A reader requesting one
        # explicitly still records its refusal through checked().
        if scope.admitted(candidate):
            yield candidate


def glob(path: Path, pattern: str) -> Iterator[Path]:
    if _scope.get() is None:
        yield from Path(path).glob(pattern)
        return
    parts = pattern.split("/")
    candidates = [Path(path)]
    for index, part in enumerate(parts):
        found = []
        for parent in candidates:
            if not is_dir(parent):
                continue
            wildcard = any(char in part for char in "*?[")
            if wildcard:
                children = (child for child in iterdir(parent)
                            if fnmatch.fnmatchcase(child.name, part))
            else:
                children = iter([parent / part])
            if index == len(parts) - 1:
                found.extend(child for child in children if wildcard or exists(child))
            else:
                found.extend(child for child in children if is_dir(child))
        candidates = found
    yield from candidates
