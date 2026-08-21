"""A campaign said "0/3 done" for hours and nothing else.

A member has no terminal. Its stage bars go to `on_step_progress`, which in
a worker goes nowhere, and the runner keeps them out of the run log on
purpose -- a file of the same line at one per cent intervals is noise. So a
member reported "Production started" and then nothing until it finished:
zero `ns/day` lines across `run.log`, `simulation.log` and
`live_events.log`, twenty minutes into production on a three-member sweep.
A healthy run and a stalled one looked identical, and the only recourse was
watching a DCD grow.

Nothing new is recorded to fix it. `live_status.json` already carries
`current_step` against `total_planned_steps`. The study's heartbeat now
reads it back.
"""
import json
from pathlib import Path

import pytest

from fastmdxplora.batch.explorer import _how_far_along


def _member(root: Path, name: str, **status):
    directory = root / name / "simulation"
    directory.mkdir(parents=True)
    if status:
        (directory / "live_status.json").write_text(json.dumps(status))
    return (name, root / name)


def test_it_reports_each_running_member(tmp_path):
    runs = [
        _member(tmp_path, "s1__random-seed-20260813",
                current_step=4_200_000, total_planned_steps=10_750_000),
        _member(tmp_path, "s1__random-seed-20260814",
                current_step=8_600_000, total_planned_steps=10_750_000),
    ]
    line = _how_far_along(runs)
    assert "39%" in line
    assert "80%" in line


def test_the_shared_prefix_is_dropped(tmp_path):
    """Every member of a sweep carries the same prefix; the axis value is
    what tells them apart, and the bar has one line to say it in."""
    runs = [_member(tmp_path, "s1__random-seed-20260813",
                    current_step=1, total_planned_steps=100)]
    line = _how_far_along(runs)
    assert "random-seed-20260813" in line
    assert "s1__" not in line


@pytest.mark.parametrize("status", [
    {},                                                  # not written yet
    {"current_step": 5},                                 # no total
    {"total_planned_steps": 100},                        # no current
    {"current_step": 5, "total_planned_steps": 0},       # nothing to divide by
    {"current_step": None, "total_planned_steps": None},  # telemetry off
])
def test_a_member_that_cannot_say_is_left_out(tmp_path, status):
    runs = [_member(tmp_path, "s1__quiet", **status)]
    assert _how_far_along(runs) == ""


def test_a_half_written_status_file_does_not_raise(tmp_path):
    """The file is rewritten while the heartbeat reads it. A heartbeat that
    fails is worse than one that says less."""
    directory = tmp_path / "s1__mid_write" / "simulation"
    directory.mkdir(parents=True)
    (directory / "live_status.json").write_text('{"current_step": 12')
    assert _how_far_along([("s1__mid_write", tmp_path / "s1__mid_write")]) == ""


def test_nothing_in_flight_adds_nothing(tmp_path):
    assert _how_far_along([]) == ""


def test_one_member_reporting_is_enough(tmp_path):
    """A campaign where one member has started and another has not still
    says what it knows."""
    runs = [
        _member(tmp_path, "s1__started",
                current_step=50, total_planned_steps=100),
        _member(tmp_path, "s1__not_yet"),
    ]
    line = _how_far_along(runs)
    assert "started 50%" in line
    assert "not_yet" not in line
