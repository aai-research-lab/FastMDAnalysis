"""A record the study pages refuse is not read, which is not missing.

A reader that falls back where a record is missing must not fall back where
a held read refused it (a link out of the study, or to a file of it that is
not one of its results). A carried-on study whose record of the join was
refused read as a study of one run: the Overview gave its first piece's
means as the production's and called them not determined. Each such
record now reads as not read, named, and no stand-in is given for it.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from tests.test_the_study_pages_read_only_the_study import (
    _a_study_joined,
    _ask,
    _link,
    _study,
    _write,
)

PUBLIC = pytest.mark.parametrize("public", [False, True], ids=["loopback", "public"])


def _energy(temperature: float, rows: int = 40) -> str:
    """A state record as OpenMM writes it, about ``temperature``."""
    return ('#"Step","Time (ps)","Potential Energy (kJ/mole)","Temperature (K)"\n'
            + "".join(f"{100 * i},{0.2 * i:.1f},{-1000.0 + i % 3},{temperature + 0.5 * (i % 4)}\n"
                      for i in range(1, rows + 1)))


def _joined(tmp_path: Path) -> Path:
    """A study carried on in two pieces of 40 samples each and joined: a
    mean over one piece is over 40, the study's over 80."""
    root = _a_study_joined(tmp_path)
    _write(root, "simulation/energy.csv", _energy(300.0))
    _write(root, "segment-001/simulation/energy.csv", _energy(300.0))
    return root


def _linked_out(root: Path, tmp_path: Path, *names: str) -> None:
    for name in names:
        path = root / name
        target = tmp_path / "outside" / name
        target.parent.mkdir(parents=True, exist_ok=True)
        path.rename(target)
        _link(path, target)


def _overview(root: Path, tmp_path: Path, *, public: bool) -> dict:
    status, body = _ask(root, tmp_path, "/api/overview", public=public)
    assert status == 200
    return json.loads(body)


def _thermodynamics_record(root: Path) -> None:
    """The thermodynamics analysis's record of the whole production, as it
    writes it (`findings.thermodynamics`)."""
    _write(root, "analysis/thermodynamics/options.json", {"findings": {"thermodynamics": {
        "temperature": {"units": "K", "mean": 300.4, "standard_error": 0.2,
                        "effective_samples": 30.0}}}})


@PUBLIC
def test_a_joined_study_s_means_are_over_both_pieces(tmp_path, public):
    """Every record read: the means are over both pieces, not one."""
    answer = _overview(_joined(tmp_path), tmp_path, public=public)
    assert answer["thermodynamics"]["means"]["temperature"]["of"] == 80
    assert "not_read" not in answer["thermodynamics"]


@PUBLIC
@pytest.mark.parametrize("refused", [
    ("joined/joined.json",),
    ("joined/joined.json", "analysis/thermodynamics/options.json"),
], ids=["join-record", "join-record-and-thermodynamics"])
def test_a_join_record_not_read_gives_no_one_piece_means(tmp_path, refused, public):
    """The record of the join refused (and the thermodynamics analysis's
    record with it): the production's means are not read, the join's
    record is named, and no piece's means stand in for the study's."""
    root = _joined(tmp_path)
    _thermodynamics_record(root)
    _linked_out(root, tmp_path, *refused)

    answer = _overview(root, tmp_path, public=public)
    thermo = answer["thermodynamics"]
    if len(refused) == 1:
        # The analysis's own record of the whole production still read.
        assert thermo["means"]["temperature"]["mean"] == pytest.approx(300.4)
        assert thermo["means"]["potential_energy"] is None
    else:
        assert thermo["means"]["temperature"] is None
    # No mean over the 40 samples of the piece read.
    assert all(mean is None or mean["of"] != 40 for mean in thermo["means"].values())
    assert thermo["not_read"] == ["joined/joined.json"]
    assert set(refused) <= set(answer["unavailable"])
    assert answer["unavailable_why"]["joined/joined.json"] == "outside"

    status, body = _ask(root, tmp_path, "/api/results", public=public)
    cards = {card["label"]: card for card in json.loads(body)["summary_cards"]}
    assert cards["Temperature"]["value"] == "Not read"
    assert cards["Temperature"]["detail"] == (
        "joined/joined.json was not read: it is a link to a file outside the study.")


@PUBLIC
@pytest.mark.parametrize("refused", [
    "joined/production.dcd", "simulation/production.dcd",
    "segment-001/simulation/energy.csv", "segment-001/simulation/checkpoint.chk.json",
])
def test_a_piece_s_record_not_read_is_named_not_missing(tmp_path, refused, public):
    """A joined trajectory, a piece's trajectory, state record or
    checkpoint's record refused: the means are not read and the record is
    named, not a study of one run, nor a piece without its record."""
    root = _joined(tmp_path)
    _linked_out(root, tmp_path, refused)

    thermo = _overview(root, tmp_path, public=public)["thermodynamics"]
    assert thermo["not_read"] == [refused]
    assert all(mean is None for mean in thermo["means"].values())


def _a_run_of_one_piece(tmp_path: Path) -> Path:
    """A study of one run whose live record sampled about 250 K and whose
    state record (`energy.csv`) about 300 K."""
    root = _study(tmp_path)
    _write(root, "simulation/live_status.json", {
        "status": "completed", "stage": "production", "current_step": 4000,
        "total_planned_steps": 4000, "last_update_timestamp": "2026-10-01T00:00:00+00:00"})
    _write(root, "simulation/live_metrics.csv", "stage,simulation_time_ns,temperature\n" + "".join(
        f"production,{0.001 * i:.3f},{250.0 + i % 3}\n" for i in range(1, 41)))
    _write(root, "simulation/energy.csv", _energy(300.0))
    return root


@PUBLIC
def test_a_state_record_not_read_is_not_stood_in_for_by_the_live_record(tmp_path, public):
    """The run's `energy.csv` refused: its means are not read, rather than
    taken from the live record's coarser samples."""
    read = _a_run_of_one_piece(tmp_path / "read")
    mean = _overview(read, tmp_path / "read", public=public)["thermodynamics"]["means"]
    assert mean["temperature"]["mean"] == pytest.approx(300.75, abs=0.5)

    root = _a_run_of_one_piece(tmp_path / "refused")
    _linked_out(root, tmp_path, "simulation/energy.csv")
    thermo = _overview(root, tmp_path / "refused", public=public)["thermodynamics"]
    assert thermo["means"]["temperature"] is None
    assert thermo["not_read"] == ["simulation/energy.csv"]


def test_outside_a_held_read_a_joined_study_is_its_pieces(tmp_path):
    """Outside a held read (the analysis phase, the command line) nothing
    is refused: a link is read as before."""
    from fastmdxplora.analysis.thermodynamics import state_record_pieces

    root = _joined(tmp_path)
    _linked_out(root, tmp_path, "joined/joined.json", "segment-001/simulation/energy.csv")
    assert len(state_record_pieces(root)) == 2


def _timed(name: str, minutes: int) -> dict:
    return {"name": name, "status": "ok", "started_at": "2026-10-01T00:00:00+00:00",
            "finished_at": f"2026-10-01T00:{minutes:02d}:00+00:00"}


@pytest.mark.parametrize("refused", ["segment-001/manifest.json", "segment-001"])
def test_a_piece_s_phases_not_read_add_no_wall_time(tmp_path, refused):
    """A carried-on study's wall time is its phases' and its pieces': with a
    piece's manifest (or the piece) refused, the others' times added are
    not the study's, and no wall time is given."""
    from fastmdxplora import record_access
    from fastmdxplora.gui.problems_view import admitted_for_problems
    from fastmdxplora.gui.simulated_time import phases_wall_seconds

    root = _a_study_joined(tmp_path)
    _write(root, "manifest.json", {"phases": [_timed("simulation", 10)]})
    _write(root, "segment-001/manifest.json", {"phases": [_timed("simulation", 20)]})
    with record_access.confined_records(root, admitted_for_problems(root)):
        assert phases_wall_seconds(root) == 30 * 60
    _linked_out(root, tmp_path, *(["segment-001/manifest.json"] if refused.endswith("json")
                                  else []))
    if refused == "segment-001":
        (root / "segment-001").rename(tmp_path / "outside-piece")
        _link(root / "segment-001", tmp_path / "outside-piece")
    with record_access.confined_records(root, admitted_for_problems(root)) as scope:
        assert phases_wall_seconds(root) is None
    assert refused in scope.denied
    # Outside a held read, nothing is refused.
    assert phases_wall_seconds(root) == 30 * 60
