"""The Report page's outline, and each number in the report traced to its record.

The report was a column of text: to find a section a person scrolled, and
a number in it (a mean, a length, a count of atoms) was taken on trust.
An outline is beside it now, its sections numbered as the PDF numbers
them, the one being read marked; and a number given with a unit, an
error, a percent or what it counts, or with a decimal point, opens the
record it was traced to: the file, the place in it, the value as
recorded and, where the report gives it in other units, as it gives it.
A number no record holds is said so, with the nearest of its kind, and
the records written again since the report are named: a report written
before its analyses were run again says numbers they no longer hold.
"""

from __future__ import annotations

import json
import os
import re
import shutil
from pathlib import Path

import pytest

from tests.test_the_frame_holds_the_study import _a_finished_study

REPORT = """# Report

## Summary

The prepared system contained 7,439 atoms and was simulated for 20 ps at 26.85 \u00b0C.

## Results

### RMSD

- Mean over the equilibrated part of the run: 0.799 \u00c5, from 8 independent samples.

### RG

- Mean over the equilibrated part of the run: 1.235 nm.
- A mean of 0.0799 nm is the RMSD's, not this section's.

## Convergence

| observable | mean |
|---|---|
| RMSD | 0.07990 |
| Rg | 1.3 |

At least 50 ps more production would determine it. The run lasted 50 ps.
"""


def _findings(mean: float, samples: float, more_ns: float | None = None) -> dict:
    found = {"unit": "nm", "mean": mean, "standard_error": float("nan"), "n_frames": 100,
             "discard": 10, "effective_samples": samples}
    if more_ns is not None:
        found["shortfall"] = {"more_ns": more_ns, "target_independent": 25.0}
    return {"options": {"selection": "protein"}, "findings": {"mean": found}}


@pytest.fixture()
def study(tmp_path) -> Path:
    root = tmp_path / "study"
    shutil.copytree(_a_finished_study(), root)
    (root / "report" / "report.md").write_text(REPORT, encoding="utf-8")
    (root / "analysis" / "rmsd" / "options.json").write_text(
        json.dumps(_findings(0.0799, 7.9, more_ns=0.05)), encoding="utf-8")
    (root / "analysis" / "rg" / "options.json").write_text(
        json.dumps(_findings(1.2345, 22.0)), encoding="utf-8")
    (root / "setup").mkdir()
    (root / "setup" / "setup_parameters.json").write_text(
        json.dumps({"n_atoms_solvated": 7439}), encoding="utf-8")
    (root / "resolved_config.yml").write_text(
        "simulation:\n  duration_ns: 0.02\n  temperature_K: 300.0\n", encoding="utf-8")
    # Everything as it was when the report was written, but the Rg
    # analysed again after it.
    written = (root / "report" / "report.md").stat().st_mtime
    for path in root.rglob("*"):
        if path.is_file():
            os.utime(path, (written - 60, written - 60))
    os.utime(root / "report" / "report.md", (written, written))
    os.utime(root / "analysis" / "rg" / "options.json", (written + 60, written + 60))
    return root


def _payload(study):
    from fastmdxplora.gui.report_page import report_payload

    said = report_payload(study)
    if said.get("rendered") != "html":
        pytest.skip("the Markdown library is not installed")
    return said


def _by_said(numbers):
    found: dict[str, list[dict]] = {}
    for number in numbers:
        found.setdefault(number["said"], []).append(number)
    return found


def test_each_number_is_traced_to_its_record_in_the_unit_the_report_gives(study):
    said = _payload(study)
    numbers = _by_said(said["numbers"])
    atoms = numbers["7,439"][0]
    assert atoms["traced"] and atoms["record"]["file"] == "setup/setup_parameters.json"
    assert atoms["record"]["key"] == "n_atoms_solvated"
    length = numbers["20 ps"][0]["record"]
    assert (length["file"], length["key"], length["unit"]) == ("resolved_config.yml",
                                                                "simulation.duration_ns", "ns")
    assert length["as_said"] == pytest.approx(20.0)
    # Kelvin recorded, Celsius said.
    assert numbers["26.85 \u00b0C"][0]["record"]["key"] == "simulation.temperature_K"
    # nm recorded, angstrom said, in its own analysis's section.
    rmsd = numbers["0.799 \u00c5"][0]
    assert rmsd["traced"] and rmsd["record"]["file"] == "analysis/rmsd/options.json"
    assert rmsd["record"]["key"] == "findings.mean.mean"
    assert numbers["8"][0]["record"]["key"] == "findings.mean.effective_samples"
    assert numbers["1.235 nm"][0]["record"]["file"] == "analysis/rg/options.json"
    # A table's row says which analysis, its column which field.
    assert numbers["0.07990"][0]["record"]["file"] == "analysis/rmsd/options.json"


def test_a_number_no_record_holds_is_said_so_with_the_nearest(study):
    said = _payload(study)
    numbers = _by_said(said["numbers"])
    # The RMSD's mean, said under the Rg: not the Rg's record, and no other
    # analysis's.
    elsewhere = numbers["0.0799 nm"][0]
    assert elsewhere["traced"] is False
    invented = numbers["1.3"][0]
    assert invented["traced"] is False
    assert invented["nearest"]["file"] == "analysis/rg/options.json"
    assert invented["nearest"]["key"] == "findings.mean.mean"
    # What a run still needs is a record of what is missing: taken for "50
    # ps more", not for the 50 ps a run lasted.
    more, lasted = numbers["50 ps"]
    assert more["traced"] and more["record"]["key"] == "findings.mean.shortfall.more_ns"
    assert lasted["traced"] is False
    assert said["records_changed"] == ["analysis/rg/options.json"]
    # The marks are around the numbers alone, and nothing else changed; the
    # report as rendered is kept for the standalone page.
    assert 'data-traced="true" role="button" tabindex="0" aria-expanded="false">7,439</span> atoms' \
        in said["marked"]
    assert "report-number" not in said["html"]
    assert re.sub(r'<span class="report-number"[^>]*>([^<]*)</span>', r"\1", said["marked"]) \
        == said["html"]


def test_a_record_linked_out_of_the_study_is_not_read_to_trace_a_number(study, tmp_path):
    """`/api/report` is read held to the study, but the numbers were traced
    with plain reads: a setup record linked out of the study gave its value
    to the report's atom count, and the outside file's time put it among
    the records changed. It is not read now, its number said as not traced
    because it was not read, and it is named as not read; also when an
    earlier trace, read without the hold, is kept."""
    from fastmdxplora.gui.report_page import report_payload
    from fastmdxplora.gui.server import _confined, _said_not_read
    from tests.test_the_study_pages_read_only_the_study import _link

    _payload(study)
    record = study / "setup" / "setup_parameters.json"
    outside = tmp_path / "outside" / "setup_parameters.json"
    outside.parent.mkdir()
    outside.write_text(json.dumps({"n_atoms_solvated": 7439}), encoding="utf-8")
    later = (study / "report" / "report.md").stat().st_mtime + 120
    os.utime(outside, (later, later))
    record.unlink()
    _link(record, outside)

    # Read without the hold first (the standalone page): it follows the link.
    plain = _by_said(report_payload(study)["numbers"])["7,439"][0]
    assert plain["traced"] and plain["record"]["file"] == "setup/setup_parameters.json"
    for _ in range(2):
        said = _said_not_read(_confined(study, report_payload))
        atoms = _by_said(said["numbers"])["7,439"][0]
        assert atoms["traced"] is False
        assert atoms["not_read"] == ["setup/setup_parameters.json"]
        assert all((n.get("record") or n.get("nearest") or {}).get("file")
                   != "setup/setup_parameters.json" for n in said["numbers"])
        assert "setup/setup_parameters.json" not in said["records_changed"]
        assert said["unavailable"] == ["setup/setup_parameters.json"]
        assert said["not_read"].startswith("setup/setup_parameters.json was not read")
        # What was read is still traced.
        assert _by_said(said["numbers"])["0.799 \u00c5"][0]["traced"]


def test_the_pdf_numbers_its_sections_as_the_outline_does():
    from fastmdxplora.report.pdf import STYLESHEET

    assert 'h2::before { content: counter(section) "  "; }' in STYLESHEET
    assert 'h3::before { content: counter(section) "." counter(subsection) "  "; }' in STYLESHEET


LOOK = """() => ({
    outline: [...document.querySelectorAll('#report-outline .outline-item')]
        .map((b) => b.querySelector('.outline-number').textContent + ' '
             + b.querySelector('.outline-words').textContent),
    numbered: [...document.querySelectorAll('#report-document h2, #report-document h3')]
        .map((h) => h.dataset.section),
    current: (document.querySelector('#report-outline .outline-item[aria-current]') || {}).textContent,
    said: document.getElementById('report-numbers-said').textContent,
    card: (document.querySelector('.report-number-card') || {}).textContent || '',
    marked: document.getElementById('report-document').dataset.marked || '',
    beside: getComputedStyle(document.getElementById('report-layout')).display,
    folded: getComputedStyle(document.querySelector('#report-outline .outline-body')).display,
})"""


def test_the_outline_is_beside_the_report_and_a_number_opens_its_record(study):
    pytest.importorskip("playwright.sync_api")
    _payload(study)
    from playwright.sync_api import sync_playwright

    from fastmdxplora.gui.server import start_dashboard_session

    session = start_dashboard_session(output=str(study), host="127.0.0.1", port=0)
    try:
        with sync_playwright() as pw:
            browser = pw.chromium.launch(args=["--enable-unsafe-swiftshader"])
            page = browser.new_page(viewport={"width": 1440, "height": 700})
            page.set_default_timeout(60000)
            errors: list[str] = []
            page.on("pageerror", lambda error: errors.append(str(error)))
            page.goto(session.url + "#report", wait_until="domcontentloaded")
            page.wait_for_selector("body.state-ready")
            page.wait_for_selector("#report-outline .outline-item", state="attached")
            # The side panel folded: the page's column is then wide enough
            # for the outline beside the report.
            page.evaluate("() => window.FastMDXFrame.setCollapsed(true)")
            page.wait_for_function("() => getComputedStyle(document.getElementById("
                                   "'report-layout')).display === 'grid'")
            first = page.evaluate(LOOK)
            page.click('#report-outline .outline-item:has-text("Convergence")')
            page.wait_for_function("() => (document.querySelector('#report-outline "
                                   ".outline-item[aria-current]') || {}).textContent === '3Convergence'")
            went = page.evaluate(LOOK)
            page.click('#report-document .report-number:text-is("0.799")')
            page.wait_for_selector(".report-number-card")
            traced = page.evaluate(LOOK)
            page.keyboard.press("Escape")
            page.click('#report-document .report-number:text-is("0.0799")')
            page.wait_for_selector('.report-number-card[data-traced="false"]')
            missing = page.evaluate(LOOK)
            page.keyboard.press("Escape")
            page.click("#report-outline .outline-mark")
            marked = page.evaluate(LOOK)
            page.set_viewport_size({"width": 1000, "height": 700})
            page.wait_for_function("() => getComputedStyle(document.getElementById("
                                   "'report-layout')).display === 'block'")
            narrow = page.evaluate(LOOK)
            page.click("#report-outline .outline-toggle")
            opened = page.evaluate(LOOK)
            browser.close()
    finally:
        session.server.shutdown()
    assert errors == []
    assert first["outline"] == ["1 Summary", "2 Results", "2.1 RMSD", "2.2 RG", "3 Convergence"]
    assert first["numbered"] == ["1", "2", "2.1", "2.2", "3"]
    assert first["beside"] == "grid" and first["folded"] == "block"
    assert first["said"].startswith("8 of 11 numbers traced to a record")
    assert "1 record changed after the report was written: analysis/rg/options.json" in first["said"]
    assert went["current"] == "3Convergence"
    assert "analysis/rmsd/options.json" in traced["card"] and "findings.mean.mean = 0.0799 nm" in traced["card"]
    assert "That is 0.799 \u00c5, as the report gives it." in traced["card"]
    assert "No record of this study holds it." in missing["card"]
    assert "analysis/rg/options.json was written again after the report" in missing["card"]
    assert marked["marked"] == "true"
    assert (narrow["beside"], narrow["folded"]) == ("block", "none")
    assert opened["folded"] == "block"


def test_the_page_says_a_number_of_a_record_not_read_is_not_traced(study, tmp_path):
    pytest.importorskip("playwright.sync_api")
    _payload(study)
    from playwright.sync_api import sync_playwright

    from fastmdxplora.gui.server import start_dashboard_session
    from tests.test_the_study_pages_read_only_the_study import _link

    record = study / "setup" / "setup_parameters.json"
    outside = tmp_path / "outside" / "setup_parameters.json"
    outside.parent.mkdir()
    outside.write_text(json.dumps({"n_atoms_solvated": 7439}), encoding="utf-8")
    record.unlink()
    _link(record, outside)
    session = start_dashboard_session(output=str(study), host="127.0.0.1", port=0)
    try:
        with sync_playwright() as pw:
            browser = pw.chromium.launch(args=["--enable-unsafe-swiftshader"])
            page = browser.new_page(viewport={"width": 1440, "height": 900})
            page.set_default_timeout(60000)
            errors: list[str] = []
            page.on("pageerror", lambda error: errors.append(str(error)))
            page.goto(session.url + "#report", wait_until="domcontentloaded")
            page.wait_for_selector("#report-outline .outline-item", state="attached")
            page.click('#report-document .report-number:text-is("7,439")')
            page.wait_for_selector('.report-number-card[data-traced="false"]')
            said = page.evaluate(LOOK)
            browser.close()
    finally:
        session.server.shutdown()
    assert errors == []
    assert said["card"].startswith("7,439×Not traced: setup/setup_parameters.json was not read, "
                                   "so whether it holds it is not known.")
    assert "No record of this study holds it." not in said["card"]
    assert "setup/setup_parameters.json was not read, so a number not traced may be one it holds." \
        in said["said"]
