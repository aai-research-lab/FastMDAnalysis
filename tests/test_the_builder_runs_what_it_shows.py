"""The Config Builder runs what it shows.

Found by the weekly test of 10-20 (area 2, New study): with a
fastmdx-defaults.yml beside the studies, the Builder's plan said 300 K
where the run used 310 K, and a value typed to set FastMDXplora's own
default back was left out of the file; `fastmdx gui` started in a folder put
its studies, read its relative paths and looked for its defaults one folder
up; Check it said "Runs." of a Config with no systems; Run as it stands went
to the Results box's old value, not the Config's own `output:`; the plan
named the pressure in atm where the run used bar; a Run without PDBFixer was
refused even for a structure already fixed, and left a folder behind; Full
config ticked still ran the short file; and the Changed card said "not set"
of a setting whose default is 2 ns. Each test here says what a person sees
or what a run is given.
"""

from __future__ import annotations

import json
import urllib.request
from pathlib import Path

import pytest

ATOM = "ATOM  {:5d}  {:<3s} {:>3s} {}{:4d}    {:8.3f}{:8.3f}{:8.3f}  1.00  0.00          {:>2s}\n"


def _peptide(n: int = 3) -> str:
    lines, serial = [], 1
    for i in range(n):
        for k, name in enumerate(("N", "CA", "C", "O")):
            lines.append(ATOM.format(serial, name, "GLY", "A", i + 1, 3.6 * i + 0.4 * k,
                                     0.0, 0.0, name[0]))
            serial += 1
    return "".join(lines) + "END\n"


@pytest.fixture(autouse=True)
def _own_folders(monkeypatch, tmp_path):
    monkeypatch.setenv("FASTMDXPLORA_CACHE_DIR", str(tmp_path / "cache"))
    monkeypatch.setenv("FASTMDXPLORA_CONFIG_DIR", str(tmp_path / "config"))
    # Your defaults are looked for up to the home folder: kept to this test.
    monkeypatch.setenv("HOME", str(tmp_path))


def _post(url: str, body: dict) -> dict:
    request = urllib.request.Request(url, data=json.dumps(body).encode(),
                                     headers={"Content-Type": "application/json"},
                                     method="POST")
    with urllib.request.urlopen(request) as response:
        return json.load(response)


def _get(url: str) -> dict:
    with urllib.request.urlopen(url) as response:
        return json.load(response)


@pytest.fixture
def home_gui(tmp_path, monkeypatch):
    """`fastmdx gui` typed in a folder of studies, with no study open,
    from a process whose own folder is elsewhere."""
    from fastmdxplora.gui import preview
    from fastmdxplora.gui.server import start_dashboard_session

    monkeypatch.setattr(preview, "_platform", lambda requested, precision: ("CPU", precision))
    studies = tmp_path / "lab" / "studies"
    studies.mkdir(parents=True)
    elsewhere = tmp_path / "elsewhere"
    elsewhere.mkdir()
    monkeypatch.chdir(elsewhere)
    session = start_dashboard_session(output=str(studies), host="127.0.0.1", port=0,
                                      home_mode=True)
    try:
        yield session, studies, f"http://127.0.0.1:{session.port}"
    finally:
        session.stop()


class TestTheFolderItWasStartedIn:

    def test_new_studies_go_in_it(self, home_gui) -> None:
        session, studies, url = home_gui
        assert Path(session.runtime.exploration_root) == studies.resolve()
        assert Path(_get(url + "/api/schema")["workspace"]) == studies.resolve()

    def test_its_defaults_are_the_form_s(self, home_gui) -> None:
        _, studies, url = home_gui
        (studies / "fastmdx-defaults.yml").write_text(
            "simulation:\n  temperature_K: 310\n", encoding="utf-8")
        schema = _get(url + "/api/schema")
        fields = {f["name"]: f for f in schema["phases"]["simulation"]["fields"]}
        assert fields["temperature_K"]["default"] == 310
        assert fields["temperature_K"].get("default_from")

    def test_a_structure_beside_it_is_previewed_as_run_will_read_it(self, home_gui) -> None:
        _, studies, url = home_gui
        (studies / "only_here.pdb").write_text(_peptide(), encoding="utf-8")
        answer = _post(url + "/api/preview-system",
                       {"system": "only_here.pdb", "include_phase": ["setup", "simulation"]})
        assert answer["ok"], answer
        assert answer["estimate"]["residues"] == 3
