"""The run's times, the Analysis page's figures and the Files page read only
the study, and say what they did not read.

The Overview's headline and The run's production came from `/api/status`,
which read the run's `simulation_parameters.json` through a link out of the
study, under the line saying that record was not read. The Analysis page's
Thermodynamics figure (`/api/figure-data`) read its record's units and
notes through a link. The Files page dropped a refused analysis's value
without a word, and a refused analysis manifest made each analysis with no
result read "done". Each is now held to the study and says what it did not
read, and why: a link out of the study, or a file of the study that is not
one of its results.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from tests.test_the_study_pages_read_only_the_study import (
    TEXT,
    _a_study_carried_on,
    _ask,
    _every_record_linked_out,
    _link,
    _study,
    _watch,
    _write,
    _OPENED,
    _WATCHED,
)

PUBLIC = pytest.mark.parametrize("public", [False, True], ids=["loopback", "public"])

#: What the outside record says the production was, which no page may say.
OUTSIDE_NS = 4242.42


def _a_run_finished(root: Path) -> Path:
    """The study's simulation finished: its live record, as the runner
    writes it, and its own record of what production ran (0.02 ns)."""
    _write(root, "simulation/live_status.json", {
        "status": "completed", "stage": "production", "current_step": 11000,
        "total_planned_steps": 11000, "timestep_fs": 2.0, "nvt_steps_planned": 500,
        "npt_steps_planned": 500, "production_steps_planned": 10000,
        "last_update_timestamp": "2026-10-01T00:00:00+00:00"})
    _write(root, "simulation/simulation_parameters.json", {
        "phase": "simulation", "duration_ns_actual": 0.02,
        "parameters": {"timestep_fs": 2.0}})
    return root


def _status(root: Path, workspace: Path, *, public: bool) -> tuple[str, dict]:
    status, body = _ask(root, workspace, "/api/status", public=public)
    assert status == 200
    return body.decode("utf-8"), json.loads(body)


@PUBLIC
def test_the_run_s_production_is_not_read_through_a_link_out(tmp_path, public):
    root = _a_run_finished(_study(tmp_path))

    # Its own record read, as before: nothing named as not read.
    _, answer = _status(root, tmp_path, public=public)
    assert answer["times"]["production_ns"] == pytest.approx(0.02)
    assert answer["times"]["equilibration_planned_ns"] == pytest.approx(0.002)
    assert answer["status"]["status"] == "completed"
    assert "unavailable" not in answer and "unavailable_why" not in answer

    outside = _write(tmp_path, "elsewhere/simulation_parameters.json", {
        "phase": "simulation", "duration_ns_actual": OUTSIDE_NS,
        "parameters": {"timestep_fs": 2.0}, "note": TEXT})
    (root / "simulation" / "simulation_parameters.json").unlink()
    _link(root / "simulation" / "simulation_parameters.json", outside)

    text, answer = _status(root, tmp_path, public=public)

    assert "4242" not in text and TEXT not in text
    assert answer["times"]["production_ns"] != pytest.approx(OUTSIDE_NS)
    assert answer["unavailable"] == ["simulation/simulation_parameters.json"]
    assert answer["unavailable_why"] == {"simulation/simulation_parameters.json": "outside"}
    # The run's state is the study's own.
    assert answer["status"]["status"] == "completed"


#: The routes beside the confined ones that a page reads the study through.
STATUS_AND_FIGURES = ("/api/status",
                      "/api/figure-data?figure=analysis/thermodynamics/thermodynamics.png")


@PUBLIC
@pytest.mark.parametrize("route", STATUS_AND_FIGURES)
def test_no_record_is_read_through_a_link_out_for_the_run_or_a_figure(tmp_path, route,
                                                                      public):
    root = _a_study_carried_on(tmp_path)
    _thermodynamics(root, units="kJ/mol")
    outside = _every_record_linked_out(root, tmp_path)
    _watch(outside)
    try:
        status, body = _ask(root, tmp_path, route, public=public)
        opened = sorted(set(_OPENED))
    finally:
        _WATCHED.clear()
    assert status == 200
    assert opened == []
    assert json.loads(body)["unavailable"]


def _thermodynamics(root: Path, *, units: str) -> None:
    """The thermodynamics analysis's table and record, as it writes them."""
    _write(root, "analysis/thermodynamics/thermodynamics.dat",
           "observable,mean,standard_error,effective_samples,standard_deviation\n"
           "potential_energy,-506500.0,25.3,30.0,140.0\n"
           "temperature,300.1,0.6,40.0,3.8\n")
    _write(root, "analysis/thermodynamics/options.json", {"findings": {"thermodynamics": {
        "potential_energy": {"units": units, "mean": -506500.0},
        "temperature": {"units": units, "mean": 300.1},
        "density": {"units": units, "not_a_measurement": TEXT}}}})


FIGURE = "/api/figure-data?figure=analysis/thermodynamics/thermodynamics.png"


@PUBLIC
def test_a_figure_s_record_linked_out_is_not_read_and_said_not_read(tmp_path, public):
    root = _study(tmp_path)
    _thermodynamics(root, units="kJ/mol")

    # Its own record read, as before.
    answer = json.loads(_ask(root, tmp_path, FIGURE, public=public)[1])
    assert answer["ok"] is True and "unavailable" not in answer
    assert answer["bars"][0]["said"][1] == "-506,500 ± 25 kJ/mol"
    assert answer["note"] == "Density: not determined"

    _thermodynamics(root, units=TEXT)
    record = root / "analysis/thermodynamics/options.json"
    outside = tmp_path / "elsewhere" / "options.json"
    outside.parent.mkdir()
    record.rename(outside)
    _link(record, outside)

    status, body = _ask(root, tmp_path, FIGURE, public=public)

    assert status == 200
    assert TEXT not in body.decode("utf-8")
    answer = json.loads(body)
    # Its bars from the study's own table, without what the record says.
    assert answer["ok"] is True
    assert [bar["label"] for bar in answer["bars"]] == ["Potential energy", "Temperature"]
    assert answer["bars"][0]["said"][1] == "-506,500 ± 25"
    assert "not determined" not in answer["note"]
    assert answer["note"] == ("analysis/thermodynamics/options.json was not read: it is a "
                              "link to a file outside the study, so what it says is not shown.")
    assert answer["unavailable"] == ["analysis/thermodynamics/options.json"]


@PUBLIC
def test_the_files_page_says_a_value_not_read(tmp_path, public):
    from tests.test_the_study_pages_read_only_the_study import OUTSIDE

    root = _study(tmp_path)
    outside = _write(tmp_path, "elsewhere/options.json", OUTSIDE["determined"])
    _link(root / "analysis/rmsd/options.json", outside)

    status, body = _ask(root, tmp_path, "/api/files-page", public=public)

    assert status == 200
    html = json.loads(body)["html"]
    assert TEXT not in html
    # Where the value was, "Not read", with why on hover; Rg's own value kept.
    assert ('<span class="files-val" data-unread title="analysis/rmsd/options.json was not '
            'read: it is a link to a file outside the study.">Not read</span>') in html
    assert '<span class="files-val">1.0300 ± 0.0040 nm</span>' in html


@PUBLIC
def test_an_analysis_manifest_not_read_leaves_no_analysis_done(tmp_path, public):
    """With its manifest refused, an analysis whose record holds no mean is
    said not read, not done: the manifest outside says it failed."""
    root = _study(tmp_path)
    manifest = {"results": {"rmsd": {"status": "failed", "message": TEXT},
                            "rg": {"status": "ok"}}}
    _write(root, "analysis/analysis_manifest.json", manifest)
    _write(root, "analysis/rmsd/options.json", {"analysis": "rmsd"})

    # Its own manifest read, each status as it says.
    rows = {row["analysis"]: row for row in json.loads(
        _ask(root, tmp_path, "/api/analysis-overview", public=public)[1])["rows"]}
    assert rows["rmsd"]["status"] == "failed" and rows["rg"]["status"] == "done"

    (root / "analysis/analysis_manifest.json").unlink()
    outside = _write(tmp_path, "elsewhere/analysis_manifest.json", manifest)
    _link(root / "analysis/analysis_manifest.json", outside)

    status, body = _ask(root, tmp_path, "/api/analysis-overview", public=public)

    assert status == 200 and TEXT not in body.decode("utf-8")
    rows = {row["analysis"]: row for row in json.loads(body)["rows"]}
    assert rows["rmsd"]["status"] == "not read"
    assert rows["rmsd"]["message"] == (
        "Its status is in analysis/analysis_manifest.json, which was not read: it is a "
        "link to a file outside the study.")
    # Rg's mean was read from its own record.
    assert rows["rg"]["status"] == "done" and rows["rg"]["mean"]["determined"] is True


@PUBLIC
def test_each_record_not_read_says_why(tmp_path, public):
    """A link out of the study is said a link out; a link to a hidden file
    of the study is said to lead to a file that is not a result."""
    from tests.test_the_study_pages_read_only_the_study import OUTSIDE

    root = _study(tmp_path)
    outside = _write(tmp_path, "elsewhere/options.json", OUTSIDE["determined"])
    _link(root / "analysis/rmsd/options.json", outside)
    hidden = root / ".kept" / "options.json"
    hidden.parent.mkdir()
    (root / "analysis/rg/options.json").rename(hidden)
    _link(root / "analysis/rg/options.json", hidden)

    answer = json.loads(_ask(root, tmp_path, "/api/overview", public=public)[1])

    assert answer["unavailable_why"] == {"analysis/rg/options.json": "leads to",
                                        "analysis/rmsd/options.json": "outside"}
    tiles = {tile["analysis"]: tile for tile in answer["tiles"]}
    assert tiles["rmsd"]["why"] == ("analysis/rmsd/options.json was not read: it is a link "
                                    "to a file outside the study.")
    assert tiles["rg"]["why"] == ("analysis/rg/options.json was not read: it leads to a file "
                                  "of the study that is not one of its results.")
