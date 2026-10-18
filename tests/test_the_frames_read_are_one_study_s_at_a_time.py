"""The frames a water map and the water sites read are kept for one study at
a time, and one request never reads them while another replaces them.

The water map asks for frames holding the map's lock and the water sites ask
holding none, so with nothing held over the frames read, one request emptied
them while another was between finding them current and reading them, and
that one failed on 'xyz' (or was handed the other study's frames).
"""

from __future__ import annotations

import threading
import time
from pathlib import Path

import numpy as np
import pytest

from fastmdxplora.gui import occupancy


def _study(root: Path, name: str, x: float) -> Path:
    import mdtraj as md

    out = root / name
    (out / "simulation").mkdir(parents=True)
    xyz = np.full((2, 3, 3), x, dtype=np.float32)  # in Å, as a DCD holds it
    with md.formats.DCDTrajectoryFile(str(out / "simulation" / "frames.dcd"), "w") as handle:
        handle.write(xyz)
    return out


class _Watched(dict):
    """The frames read, slowed where the race lies: a request that finds the
    frames current waits there for another to begin replacing them, and the
    one replacing them pauses between emptying and filling."""

    def __init__(self, checking: threading.Event, emptied: threading.Event,
                 current_for: str) -> None:
        super().__init__()
        self.checking, self.emptied, self.current_for = checking, emptied, current_for

    def get(self, key, default=None):
        value = super().get(key, default)
        if (key == "key" and value == self.current_for
                and threading.current_thread().name == "first"):
            self.checking.set()
            self.emptied.wait(3.0)
        return value

    def clear(self) -> None:
        super().clear()
        self.emptied.set()
        time.sleep(0.3)


def test_a_study_s_frames_are_not_taken_while_another_s_are_read(tmp_path,
                                                                   monkeypatch) -> None:
    first = _study(tmp_path, "first", 10.0)
    second = _study(tmp_path, "second", 20.0)
    path = first / "simulation" / "frames.dcd"
    current = f"{path}:{path.stat().st_mtime_ns}"

    checking, emptied = threading.Event(), threading.Event()
    watched = _Watched(checking, emptied, current)
    monkeypatch.setattr(occupancy, "_FRAMES_READ", watched)
    occupancy._FRAMES_READ.update(key=current,
                                  xyz=np.full((2, 3, 3), 1.0, dtype=np.float32))

    results: dict[str, object] = {}

    def read(name: str, out: Path) -> None:
        try:
            results[name] = occupancy._frame_xyz(out, 0)
        except Exception as exc:  # noqa: BLE001 - the failure is what is tested
            results[name] = exc

    reader = threading.Thread(target=read, args=("first", first), name="first")
    reader.start()
    assert checking.wait(5.0)
    replacer = threading.Thread(target=read, args=("second", second), name="second")
    replacer.start()
    reader.join(10.0)
    replacer.join(10.0)

    assert not isinstance(results["first"], Exception), results["first"]
    assert not isinstance(results["second"], Exception), results["second"]
    assert np.allclose(results["first"], 1.0)
    assert np.allclose(results["second"], 2.0)


@pytest.fixture(autouse=True)
def _frames_read_left_as_found():
    kept = dict(occupancy._FRAMES_READ)
    yield
    occupancy._FRAMES_READ.clear()
    occupancy._FRAMES_READ.update(kept)
