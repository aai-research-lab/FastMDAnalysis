"""A mean noted for the solute's own periodic image says so on every page
that shows it, and to the Agent, beside its verdict, which it does not
change.

The report and the Report page said of an RMSD "may be biased: the solute
met its own periodic image in 108 of 108 frames used", while the Overview
tile and the Analysis page's table said "Determined 0.1046 ± 0.0013 nm" and
nothing else, and the Agent said "determined" with no caveat: one mean, two
verdicts. Every reader now reads the note from the record through one
reader (`fastmdxplora.image_note`).
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from fastmdxplora.image_note import said_short
from tests.test_the_analysis_page_reads_as_a_whole import _analysis, _ar1, _manifest

N = 2000
#: Frames 0, 5, 10, ... came within the cutoff of the solute's own image.
EVERY = 5


def _finding() -> tuple[dict, np.ndarray]:
    """The image check's finding and series: within 1 nm of its image in
    every fifth frame, 0.42 nm the closest."""
    series = np.full(N, 2.0)
    series[::EVERY] = 0.8
    series[EVERY * 7] = 0.42
    finding = {
        "cutoff_nm": 1.0, "cutoff_from": "taken because no study was given", "close_nm": 0.5,
        "frames": N, "nearest_nm": 0.42, "nearest_frame": EVERY * 7,
        "frames_under_cutoff": N // EVERY, "frames_under_close": 1,
        "reaches_its_image": True, "solute": "the macromolecules (10 residues in 1 chain); "
                                            "water, ions and lipids left out",
        "said": ("The solute came within 0.420 nm of its own periodic image, inside the 1 nm "
                 f"cutoff in {N // EVERY} of {N} frames (20%) and in contact with it (under "
                 "0.5 nm) in 1 (0.05%). The closest was at frame 35, between ALA1 CA and "
                 "the image of ALA10 CA (heavy atoms only)."),
    }
    return finding, series


def _noted(root: Path, name: str, finding: dict, series: np.ndarray) -> dict:
    """The analysis's record noted as the analysis code notes it (one
    function, `note_on_frames`), over the frames its mean uses."""
    from fastmdxplora.analysis.image_check import note_numbers, note_on_frames

    path = root / "analysis" / name / "options.json"
    record = json.loads(path.read_text())
    findings = record["findings"]
    for key, mean in list(findings.items()):
        if isinstance(mean, dict) and "n_frames" in mean:
            note = note_on_frames(finding, series, np.arange(mean["discard"], N),
                                  analysis=name)
            mean["solute_image"] = note_numbers(note)
            if key == "mean":
                findings["solute_image"] = note
    path.write_text(json.dumps(record))
    return findings


@pytest.fixture
def study(tmp_path: Path) -> Path:
    root = tmp_path / "study"
    _manifest(root, {"image_check": {"status": "ok"}, "rmsd": {"status": "ok"},
                     "sasa": {"status": "ok"}})
    finding, series = _finding()
    folder = root / "analysis" / "image_check"
    folder.mkdir(parents=True)
    (folder / "image_check.dat").write_text("".join(f"{v}\n" for v in series))
    (folder / "options.json").write_text(json.dumps(
        {"analysis": "image_check", "selection": "protein",
         "findings": {"solute_image": finding}}))
    _analysis(root, "rmsd", 0.30 + 0.002 * _ar1(N, 0.5, 1))
    _analysis(root, "sasa", 91.0 + 0.3 * _ar1(N, 0.5, 3), unit="nm²",
              others={"hydrophobic_sasa": 45.0 + 0.1 * _ar1(N, 0.5, 4)})
    for name in ("rmsd", "sasa"):
        _noted(root, name, finding, series)
    return root


def _expected(root: Path, name: str, key: str = "mean") -> str:
    record = json.loads((root / "analysis" / name / "options.json").read_text())
    return said_short(record["findings"][key]["solute_image"])


class TestTheShortForm:
    def test_it_says_the_frames_the_cutoff_and_the_solute_s_own_image(self, study):
        said = _expected(study, "rmsd")
        assert said.startswith("may be biased: ")
        assert "frames used" in said and "came within 1 nm of the solute's own image" in said


class TestTheAnalysisPage:
    def test_a_noted_mean_keeps_its_verdict_and_says_the_note(self, study):
        from fastmdxplora.gui.analysis_overview import overview_of

        rows = {row["analysis"]: row for row in overview_of(study)["rows"]}
        rmsd = rows["rmsd"]["mean"]
        assert rmsd["determined"] is True and rmsd["why"] is None
        assert rmsd["noted"] == _expected(study, "rmsd")
        assert rmsd["noted_why"].startswith("The solute came within 1 nm")
        hydrophobic = next(q for q in rows["sasa"]["quantities"]
                           if q["key"] == "hydrophobic_sasa")
        assert hydrophobic["noted"] == _expected(study, "sasa", "hydrophobic_sasa")

    def test_the_image_check_s_row_carries_its_finding_as_a_sentence(self, study):
        from fastmdxplora.gui.analysis_overview import overview_of

        row = next(r for r in overview_of(study)["rows"] if r["analysis"] == "image_check")
        assert row["finding"] == (f"The solute came within 0.420 nm of its own periodic image, "
                                  f"inside the 1 nm cutoff in {N // EVERY} of {N} frames (20%) "
                                  "and in contact with it (under 0.5 nm) in 1 (0.05%).")
        assert row["finding_why"].startswith(row["finding"])
        assert row["atoms_said"] == "heavy atoms of the macromolecules (10 residues in 1 chain)"
        assert row["quantities"] == []

    def test_the_convergence_says_it_beside_the_recorded_mean(self, study):
        from fastmdxplora.gui.analysis_overview import convergence_payload

        recorded = convergence_payload(study, "rmsd")["recorded"]
        assert recorded["noted"] == _expected(study, "rmsd")

    def test_the_series_carries_it_with_its_mean(self, study):
        from fastmdxplora.gui.series import series_payload

        assert series_payload(study, "rmsd")["mean"]["noted"] == _expected(study, "rmsd")


class TestTheOverview:
    def test_a_noted_tile_says_it_and_stays_determined(self, study):
        from fastmdxplora.gui.overview_view import overview_payload

        tiles = {t["label"]: t for t in overview_payload(study)["tiles"]}
        rmsd = next(t for t in tiles.values() if t["analysis"] == "rmsd")
        assert rmsd["determined"] is True
        assert rmsd["noted"] == _expected(study, "rmsd")
        assert rmsd["noted_why"].startswith("The solute came within 1 nm")


class TestTheCardsAndTheReportPage:
    def test_the_card_caption_is_the_same_note(self, study):
        from fastmdxplora.gui.report_dashboard import _what_the_analysis_found

        caption = _what_the_analysis_found(study / "analysis" / "rmsd" / "rmsd.dat")
        assert caption.endswith("; " + _expected(study, "rmsd"))

    def test_the_image_check_s_card_is_titled_and_captioned_with_its_finding(self, study):
        from fastmdxplora.gui.report_dashboard import (
            _dashboard_summaries,
            _figure_title_from_path,
        )

        assert _figure_title_from_path(study / "analysis" / "image_check" /
                                       "image_check.png") == "Image check"
        summary = _dashboard_summaries(study)["Image check"].summary
        assert summary.startswith("The solute came within 0.420 nm") and summary.endswith("%).")

    def test_the_report_page_s_table_says_it_beside_the_verdict(self, study):
        from fastmdxplora.gui.report_dashboard import _metric_rows, _render_metric_row

        row = next(r for r in _metric_rows(study, {}) if r.metric == "RMSD")
        assert row.status == "Determined" and row.noted == _expected(study, "rmsd")
        from html import escape

        assert escape(_expected(study, "rmsd")) in _render_metric_row(row)


class TestTheAgent:
    def test_its_reading_of_the_analyses_says_it_and_the_paragraph_once(self, study):
        from fastmdxplora.gui.agent_panel import _results_summary

        said = _results_summary(study)
        rmsd = next(line for line in said.splitlines() if line.strip().startswith("rmsd:"))
        assert _expected(study, "rmsd") in rmsd
        assert "the image note: The solute came within 1 nm" in rmsd
        asked = _results_summary(study, first="rmsd")
        line = next(line for line in asked.splitlines() if line.strip().startswith("rmsd:"))
        # The analysis asked about keeps the paragraph whole on its line.
        assert "the image note: The solute came within 1 nm" in line
        paragraph = line.split("the image note: ")[1]
        assert "; the closest, " in paragraph and "[" not in paragraph

    def test_its_answer_from_the_records_says_it(self, study):
        from fastmdxplora.gui.records_answer import _means

        line = next(line for line in _means(study) if line.startswith("RMSD:"))
        assert _expected(study, "rmsd") in line

    def test_a_study_card_s_mean_says_it(self, study):
        from fastmdxplora.gui.workspace import _means
        from fastmdxplora.workspace_studies import mean_said

        side = dict(_means(study))["rmsd"]
        assert side["noted"] == _expected(study, "rmsd")
        assert mean_said(side).endswith("; " + _expected(study, "rmsd"))


# --------------------------------------------------------------------------
# In a browser: the Overview tile, the Analysis page's table and the cell
# over its card give the same note for one mean.
# --------------------------------------------------------------------------

def test_one_note_on_the_overview_the_table_and_the_card(tmp_path) -> None:
    pytest.importorskip("playwright.sync_api")
    from playwright.sync_api import sync_playwright

    from fastmdxplora.gui.server import start_dashboard_session
    from tests.test_the_drawing_scripts_run_in_a_browser import _write_study
    from tests.test_the_overview_leads_with_what_was_determined import (
        _metrics,
        _phases,
        _status,
    )

    root = tmp_path / "study"
    _write_study(root)
    (root / "simulation" / "energy.csv").unlink(missing_ok=True)
    _manifest(root, {"image_check": {"status": "ok"}, "rmsd": {"status": "ok"}})
    finding, series = _finding()
    folder = root / "analysis" / "image_check"
    folder.mkdir(parents=True)
    (folder / "image_check.dat").write_text("".join(f"{v}\n" for v in series))
    (folder / "options.json").write_text(json.dumps(
        {"analysis": "image_check", "selection": "protein",
         "findings": {"solute_image": finding}}))
    _analysis(root, "rmsd", 0.30 + 0.002 * _ar1(N, 0.5, 1))
    _noted(root, "rmsd", finding, series)
    # Each analysis's figure, which its card on the Analysis page shows.
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    for name in ("rmsd", "image_check"):
        figure = plt.figure(figsize=(2, 1.5))
        plt.plot([0, 1], [0, 1])
        figure.savefig(root / "analysis" / name / f"{name}.png")
        plt.close(figure)
    _metrics(root)
    _status(root)
    _phases(root)
    expected = _expected(root, "rmsd")

    session = start_dashboard_session(output=str(root), host="127.0.0.1", port=0)
    try:
        with sync_playwright() as pw:
            browser = pw.chromium.launch(args=["--enable-unsafe-swiftshader"])
            found = {}
            for width, scheme in ((1440, "light"), (390, "dark")):
                page = browser.new_page(viewport={"width": width, "height": 900},
                                        color_scheme=scheme)
                page.set_default_timeout(60000)
                errors: list[str] = []
                page.on("pageerror", lambda error, seen=errors: seen.append(str(error)))
                page.goto(session.url + "#overview", wait_until="domcontentloaded")
                page.wait_for_selector('#overview-tiles [data-tile-analysis="rmsd"] .tile-noted')
                tile = page.text_content('#overview-tiles [data-tile-analysis="rmsd"] .tile-noted')
                badge = page.text_content('#overview-tiles [data-tile-analysis="rmsd"] .tile-badge')
                page.goto(session.url + "#analysis", wait_until="domcontentloaded")
                page.wait_for_selector('#analysis-results-table tr[data-analysis="rmsd"] '
                                       '.analysis-results-noted summary')
                table = page.text_content('#analysis-results-table tr[data-analysis="rmsd"] '
                                          '.analysis-results-noted summary')
                status = page.text_content('#analysis-results-table tr[data-analysis="rmsd"] '
                                           '.analysis-results-status')
                check = page.text_content('#analysis-results-table '
                                          'tr[data-analysis="image_check"] td')
                page.wait_for_selector(".cell-noted")
                cell = page.text_content(".cell-noted")
                wide = page.evaluate("() => document.documentElement.scrollWidth "
                                     "- document.documentElement.clientWidth")
                found[(width, scheme)] = (tile, badge, table, status, check, cell, wide, errors)
                page.close()
            browser.close()
    finally:
        session.server.shutdown()
    for (width, scheme), (tile, badge, table, status, check, cell, wide,
                          errors) in found.items():
        assert not errors, (width, scheme, errors)
        assert tile == table == cell == expected, (width, scheme, tile, table, cell)
        assert badge == "Determined" and status.strip() == "Determined"
        assert check.startswith("The solute came within 0.420 nm"), check
        assert wide <= 0, (width, scheme, wide)
