"""The Overview and Analysis pages read only the study's own records.

A study folder can hold a link that leads out of it (a study shared as a
download, say). What the link leads to is never read for a page, on this
computer or on a public bind: the record is said to be unavailable instead,
and the pages keep what the study's own records say.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from fastmdxplora.statistics import MEAN_RULES

#: Frames the study's analyses read.
N = 40
#: What the records outside the study say, which no page may show.
TEXT = "OUTSIDE-TEXT-not-for-any-page"
VALUE = -4242.42
#: How a page could say the value: as recorded, to four figures, with its error.
VALUE_SAID = ("-4242.4", "-4,242", "-4242 ")

#: A record outside the study, in the form the analyses write
#: (`findings.mean`, `statistics.mean_record`): one a mean was determined
#: of, and one whose mean was withheld for a reason.
OUTSIDE = {
    "determined": {"findings": {"mean": {
        "mean": VALUE, "standard_error": 4.2, "standard_deviation": 30.0,
        "effective_samples": 30.0, "discard": 4, "n_frames": N, "rules": MEAN_RULES,
        "unit": TEXT}}},
    "withheld": {"findings": {"mean": {
        "mean": VALUE, "standard_deviation": 30.0, "effective_samples": 3.0,
        "discard": 4, "n_frames": N, "rules": MEAN_RULES, "unit": TEXT,
        "not_a_measurement": TEXT}}},
}

#: The routes the Overview and the Analysis page read the study through.
ROUTES = ("/api/overview", "/api/results", "/api/series", "/api/series?analysis=rmsd",
          "/api/analysis-overview", "/api/convergence?analysis=rmsd", "/api/files-page")


def _write(root: Path, relative: str, value: object) -> Path:
    path = root / relative
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(value if isinstance(value, str) else json.dumps(value), encoding="utf-8")
    return path


def _link(path: Path, target: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    try:
        path.symlink_to(target)
    except OSError:
        pytest.skip("this system cannot create symbolic links")


def _series(start: float) -> str:
    return "# whitespace-delimited\n" + "\n".join(
        f"{start + 0.01 * (i % 7):.6f}" for i in range(N)) + "\n"


def _study(tmp_path: Path) -> Path:
    """A study analysed for RMSD and Rg, each with its data file, Rg with
    its record, recorded on a 10 ps clock."""
    root = tmp_path / "study"
    _write(root, "manifest.json", {"phases": [{"name": "analysis", "status": "ok"}]})
    _write(root, "analysis/analysis_manifest.json", {
        "phase": "analysis", "n_frames": N,
        "load_kwargs": {"stride": None, "first": None, "saving_interval_ps": 10.0},
        "results": {"rmsd": {"status": "ok"}, "rg": {"status": "ok"}}})
    _write(root, "analysis/rmsd/rmsd.dat", _series(0.1))
    _write(root, "analysis/rg/rg.dat", _series(1.0))
    _write(root, "analysis/rg/options.json", {"analysis": "rg", "findings": {"mean": {
        "mean": 1.03, "standard_error": 0.004, "standard_deviation": 0.02,
        "effective_samples": 30.0, "discard": 4, "n_frames": N, "rules": MEAN_RULES,
        "unit": "nm"}}})
    return root


def _ask(root: Path, workspace: Path, route: str, *, public: bool) -> tuple[int, bytes]:
    """``route`` as a page asks it, on loopback or as a public bind
    (`--host 0.0.0.0`) answers it."""
    from fastmdxplora.gui.exploration import DashboardRuntime
    from tests.test_a_public_dashboard_answers_only_what_it_lists import _ask as ask
    from tests.test_a_public_dashboard_answers_only_what_it_lists import _serving

    runtime = DashboardRuntime(workspace_root=root, exploration_root=workspace,
                               active_root=root)
    with _serving(runtime, allow_control=not public) as url:
        return ask(url + route)


def _shows_nothing_from_outside(body: bytes) -> None:
    text = body.decode("utf-8")
    assert TEXT not in text
    for said in VALUE_SAID:
        assert said not in text, said

    def numbers(value: object):
        if isinstance(value, dict):
            for inner in value.values():
                yield from numbers(inner)
        elif isinstance(value, list):
            for inner in value:
                yield from numbers(inner)
        elif isinstance(value, (int, float)) and not isinstance(value, bool):
            yield value

    assert all(abs(number - VALUE) > 1e-6 for number in numbers(json.loads(text)))


REFUSED = "analysis/rmsd/options.json was not read"


@pytest.mark.parametrize("public", [False, True], ids=["loopback", "public"])
@pytest.mark.parametrize("route,kind", [
    (route, kind) for route in ROUTES for kind in sorted(OUTSIDE)
    # The Files page gives a value only where it was determined.
    if not (route == "/api/files-page" and kind == "withheld")])
def test_a_record_linked_out_of_the_study_reaches_no_page(tmp_path, route, kind, public):
    root = _study(tmp_path)
    outside = _write(tmp_path, "elsewhere/options.json", OUTSIDE[kind])
    _link(root / "analysis/rmsd/options.json", outside)

    status, body = _ask(root, tmp_path, route, public=public)

    assert status == 200
    _shows_nothing_from_outside(body)
    answer = json.loads(body)
    # Said as not read, where the route says anything of the record.
    if route == "/api/overview":
        tiles = {tile["analysis"]: tile for tile in answer["tiles"]}
        assert tiles["rmsd"]["said"] == "record unavailable"
        assert REFUSED in tiles["rmsd"]["why"]
        assert tiles["rg"]["determined"] is True
        assert answer["unavailable"]
    elif route == "/api/results":
        rows = {row["metric"]: row for row in answer["metric_rows"]}
        assert rows["RMSD"]["average"] == "record unavailable"
        assert REFUSED in rows["RMSD"]["why"]
        assert answer["unavailable"]
    elif route == "/api/series":
        assert answer["unavailable"]
    elif route == "/api/analysis-overview":
        rows = {row["analysis"]: row for row in answer["rows"]}
        rmsd = rows["rmsd"]["quantities"][0]
        assert rows["rmsd"]["kind"] == "mean"
        assert rmsd["said"] == "record unavailable" and rmsd["determined"] is False
        assert REFUSED in rmsd["why"]
        assert rows["rg"]["mean"]["determined"] is True
    elif route == "/api/files-page":
        # The study's own value is still given beside its analysis.
        assert "1.0300 ± 0.0040 nm" in answer["html"]
    else:
        assert answer["ok"] is False
        assert REFUSED in answer["reason"]
        assert answer["unavailable"] == ["analysis/rmsd/options.json"]


#: The run's records the time of each analysed frame is read from: the
#: analysis manifest's clock and the loader's own times, and the run's
#: record and config, which say whether the integrator chose its own step.
CLOCKS = {
    "analysis_manifest": ("analysis/analysis_manifest.json", {
        "results": {"rmsd": {"status": "ok"}, "rg": {"status": "ok"}},
        "load_kwargs": {"saving_interval_ps": 1.0},
        "frame_times_ps": [987654.0 + i for i in range(N)]}),
    "simulation_parameters": ("simulation/simulation_parameters.json", {
        "parameters": {"integrator": "variable_langevin"}}),
    "resolved_config": ("resolved_config.yml", "simulation:\n  integrator: variable_langevin\n"),
}


@pytest.mark.parametrize("public", [False, True], ids=["loopback", "public"])
@pytest.mark.parametrize("record", sorted(CLOCKS))
@pytest.mark.parametrize("route", ["/api/analysis-overview", "/api/overview",
                                   "/api/series?analysis=rg"])
def test_a_clock_linked_out_of_the_study_times_no_frame(tmp_path, route, record, public):
    """The frames' times come from the study's own records: a record linked
    out of it neither gives them nor takes them away."""
    root = _study(tmp_path)
    relative, content = CLOCKS[record]
    if record == "analysis_manifest":
        (root / relative).unlink()
    outside = _write(tmp_path, "elsewhere/" + Path(relative).name, content)
    _link(root / relative, outside)

    status, body = _ask(root, tmp_path, route, public=public)

    assert status == 200
    text = body.decode("utf-8")
    assert "987654" not in text and "987.65" not in text
    assert "variable_langevin" not in text
    answer = json.loads(body)
    if route == "/api/series?analysis=rg":
        times, start = answer["x"], answer["mean"]["from_x"]
        said = answer["x_label"]
    else:
        rows = answer["rows"] if route == "/api/analysis-overview" else answer["tiles"]
        rg = next(row for row in rows if row["analysis"] == "rg")
        mean = rg["mean"] if route == "/api/analysis-overview" else rg
        start, times, said = mean["from_ns"], None, None
    if record == "analysis_manifest":
        # Its own manifest refused, the study records no clock: said by
        # frame, not by time.
        assert start in (None, 4.0)
        if said is not None:
            assert said == "Frame" and times[0] == 0.0
    else:
        # Frame 4 of a trajectory saved every 10 ps, as the study's own
        # manifest says, the outside record not read.
        assert start == pytest.approx(0.05)
        if said is not None:
            assert said == "Time (ns)" and times[0] == pytest.approx(0.01)
