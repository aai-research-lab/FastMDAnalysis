"""A piece stopped before its first checkpoint adds nothing to what was done.

What a study had done was the sum of its pieces, and a piece with no
checkpoint was counted at the length it was asked for: an extension of
0.2 ns killed before its checkpoint added 0.2 ns, though its frames are
cut away and the span is run again from the checkpoint before it.
"""

from __future__ import annotations

import os
import subprocess
import sys

import pytest

from fastmdxplora.simulation.resume import production_done_ns
from tests.test_a_study_is_extended_in_place import _study
from tests.test_an_extended_study_goes_on_from_its_last_segment import _segment


def test_a_piece_with_no_checkpoint_counts_nothing() -> None:
    root = _study()
    _segment(root, 1)
    killed = _segment(root, 2, duration=0.2, finished=False)
    for written in (killed / "simulation").glob("checkpoint*"):
        written.unlink()
    assert production_done_ns(root) == pytest.approx(0.6)


def test_a_piece_with_its_checkpoint_counts_what_it_ran() -> None:
    root = _study()
    _segment(root, 1)
    _segment(root, 2, duration=0.2, done_steps=25_000, finished=False)
    assert production_done_ns(root) == pytest.approx(0.65)


def test_a_piece_killed_as_its_checkpoint_was_written_counts_nothing() -> None:
    """Killed between writing the checkpoint and its record of itself, a
    piece was counted at the whole length it was asked for."""
    from fastmdxplora.simulation.runner import CHECKPOINT_SIDECAR_SUFFIX

    root = _study()
    _segment(root, 1)
    killed = _segment(root, 2, duration=0.2, done_steps=25_000, finished=False)
    checkpoint = killed / "simulation" / "checkpoint.chk"
    checkpoint.with_suffix(checkpoint.suffix + CHECKPOINT_SIDECAR_SUFFIX).unlink()
    assert checkpoint.is_file()
    assert production_done_ns(root) == pytest.approx(0.6)


@pytest.mark.skipif(os.name == "nt", reason="the lock is held here with fcntl")
def test_a_study_carried_on_elsewhere_is_not_carried_on_twice(monkeypatch) -> None:
    """A second resume, or an extension, of a study another process is
    carrying on was not refused until that process had written its record:
    it would have run the same stretch into the same next piece."""
    from fastmdxplora import FastMDXplora
    from fastmdxplora.simulation.resume import CARRYING_ON_FILE, extend_study, resume_study

    root = _study(finished=False)
    held = subprocess.Popen(
        [sys.executable, "-c",
         "import fcntl, os, sys, time\n"
         "handle = os.open(sys.argv[1], os.O_RDWR | os.O_CREAT)\n"
         "fcntl.flock(handle, fcntl.LOCK_EX)\n"
         "print('held', flush=True)\n"
         "time.sleep(120)\n", str(root / CARRYING_ON_FILE)],
        stdout=subprocess.PIPE, text=True)
    ran = []
    monkeypatch.setattr(FastMDXplora, "explore", lambda self, *a, **k: ran.append(self) or [])
    try:
        assert held.stdout.readline().strip() == "held"
        answer = resume_study(root)
        extended = extend_study(root, more_ns=0.0002)
    finally:
        held.kill()
        held.wait()
    assert answer == {"study": str(root.resolve()), "ok": False, "did": "nothing",
                      "error": "A run of this study is still going. Resume it once it has "
                               "stopped."}
    assert extended["ok"] is False and "still going" in extended["error"], extended
    assert ran == []


def test_its_fix_keeps_the_study_s_production() -> None:
    """Its fix said it stopped before production and was to be run again,
    though the study's production up to the checkpoint before it is kept
    and carried on from."""
    from fastmdxplora.remedies import STOPPED_REFUSAL, _resume

    root = _study()
    _segment(root, 1)
    killed = _segment(root, 2, duration=0.2, finished=False)
    for written in (killed / "simulation").glob("checkpoint*"):
        written.unlink()
    fix = _resume(STOPPED_REFUSAL, "The run was stopped.", "the study", [killed], root).fix
    assert "before production" not in fix, fix
    assert "Production already written is kept" in fix, fix


def test_production_off_the_frame_grid_is_not_run_again_from_its_start(monkeypatch) -> None:
    """A study of one piece whose checkpoint is off its frames' grid was read
    as nothing done, and `fastmdx resume` ran it again from its start over
    the production it had written."""
    from fastmdxplora import FastMDXplora
    from fastmdxplora.simulation.resume import continuation_of, resume_study
    from fastmdxplora.simulation.runner import write_checkpoint_sidecar

    root = _study(finished=False)
    write_checkpoint_sidecar(root / "simulation" / "checkpoint.chk", stage="production",
                             step=100_100, ensemble="npt", temperature_K=300.0,
                             timestep_fs=2.0, study=str(root), finished=False,
                             trajectory_interval_steps=500)
    ran = []
    monkeypatch.setattr(FastMDXplora, "explore", lambda self, *a, **k: ran.append(self) or [])
    plan = continuation_of(root)
    answer = resume_study(root)
    assert not plan.possible and plan.production_done_ns == pytest.approx(0.2002)
    assert answer["ok"] is False and ran == [], answer
    assert "0.200 ns of production is written" in answer["error"], answer


def test_the_overview_s_speed_is_the_last_piece_run(monkeypatch, tmp_path) -> None:
    """The Overview's production speed was read from the last `segment-*`
    folder by name: a piece set aside, or `segment-9` beside `segment-010`,
    was read in place of the piece run last."""
    from fastmdxplora import remedies
    from fastmdxplora.gui.overview_view import _production_speed

    root = tmp_path / "study"
    for name in ("segment-001", "segment-001-stopped-before-its-checkpoint",
                 "segment-9", "segment-010", "segment-010-stopped-before-its-checkpoint"):
        (root / name / "simulation").mkdir(parents=True)
        (root / name / "simulation" / "live_status.json").write_text("{}")
    read = []
    monkeypatch.setattr(remedies, "_speed",
                        lambda folder, study: (read.append(folder.name) or 86400.0, "CPU"))
    assert _production_speed(root) == pytest.approx(1.0)
    assert read == ["segment-010"]


def test_a_folder_that_is_not_a_study_is_left_as_it_was(tmp_path) -> None:
    """Resumed or extended, a folder that is not a study was said to be so,
    with the study's lock file left in it."""
    from fastmdxplora.simulation.resume import extend_study, resume_study

    answer = resume_study(tmp_path)
    extended = extend_study(tmp_path, more_ns=0.0002)
    assert answer["ok"] is False and "is not a study" in answer["error"], answer
    assert extended["ok"] is False, extended
    assert list(tmp_path.iterdir()) == []


def test_the_lock_held_through_the_report_is_left_out_of_its_bundle(tmp_path) -> None:
    """The report's bundle read the lock file held as the report is made,
    which Windows does not let any handle read."""
    from fastmdxplora.report.bundle import _iter_project_files
    from fastmdxplora.simulation.resume import CARRYING_ON_FILE

    (tmp_path / CARRYING_ON_FILE).write_text("{}", encoding="utf-8")
    (tmp_path / "resolved_config.yml").write_text("{}", encoding="utf-8")
    files = _iter_project_files(tmp_path, tmp_path / "project_bundle.zip")
    assert [path.name for path in files] == ["resolved_config.yml"]
