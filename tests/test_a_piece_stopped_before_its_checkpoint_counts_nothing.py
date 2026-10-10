"""A piece stopped before its first checkpoint adds nothing to what was done.

What a study had done was the sum of its pieces, and a piece with no
checkpoint was counted at the length it was asked for: an extension of
0.2 ns killed before its checkpoint added 0.2 ns, though its frames are
cut away and the span is run again from the checkpoint before it.
"""

from __future__ import annotations

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
