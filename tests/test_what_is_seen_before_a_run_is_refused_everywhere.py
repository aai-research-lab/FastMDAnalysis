"""What can be found wrong before anything is built is refused by every door.

Found by the review of the weekly test's fixes, 10-20: a system setup could
not read (a structure file that is not there, a name that is neither a file
nor a PDB ID, a sequence) and a selection MDTraj could not parse were
refused by `fastmdx explore --dry-run` alone. The GUI's Check it said
"Runs.", an AI app's `check_study` "Accepted by the validator", and a run
without `--dry-run` went ahead and left its folder behind. And the dry run
refused two valid studies: one simulating a system prepared before
(`setup_from`), whose structure file setup never reads, and one giving
`water_sites` its own default, `site_selection: auto`.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest
import yaml

sys.path.insert(0, str(Path(__file__).parent))
from test_report_wiring import _make_traj_files  # noqa: E402

from fastmdxplora.cli.main import main  # noqa: E402

BAD = "protein and nam CA"

#: Each: the systems block, what the refusal says, and its code.
SEEN = [
    ([{"id": "a", "system": "not_here.pdb"}], "No structure at 'not_here.pdb'",
     "environment.path.not_found"),
    ([{"id": "a", "system": "1L2YA"}], "Could not classify system input '1L2YA'",
     "setup.input.unrecognised"),
    ([{"id": "a", "system": "MKTAYIAKQRQISFVKSHFSRQ"}], "reads as an amino-acid sequence",
     "environment.backend.missing"),
]
SELECTION = ("is not an MDTraj selection", "config.option.not_permitted")


@pytest.fixture(autouse=True)
def _own_folders(monkeypatch, tmp_path):
    monkeypatch.setenv("FASTMDXPLORA_CACHE_DIR", str(tmp_path / "cache"))
    monkeypatch.setenv("FASTMDXPLORA_CONFIG_DIR", str(tmp_path / "settings"))
    monkeypatch.setenv("HOME", str(tmp_path))


def _write(path: Path, data: dict) -> Path:
    path.write_text(yaml.safe_dump(data), encoding="utf-8")
    return path


def _cases():
    for systems, says, code in SEEN:
        yield {"systems": systems}, says, code
    yield ({"systems": [{"id": "a", "system": "1L2Y"}], "analysis": {"select_atoms": BAD}},
           *SELECTION)
    yield ({"systems": [{"id": "a", "system": "1L2Y",
                         "analysis": {"options": {"rmsd": {"selection": BAD}}}}]},
           *SELECTION)


CASES = list(_cases())
IDS = ["missing-file", "not-an-id", "sequence", "selection", "a-system-s-selection"]


@pytest.mark.parametrize("config, says, code", CASES, ids=IDS)
def test_a_run_is_refused_before_it_writes_anything(tmp_path, monkeypatch, capsys,
                                                     config, says, code):
    monkeypatch.chdir(tmp_path)
    # Setup alone, from a structure here, so a check that let it through
    # fails at once rather than fetching and simulating.
    _make_traj_files(tmp_path / "made")
    here = str(tmp_path / "made" / "simulation" / "topology.pdb")
    config = {**config, "systems": [
        {**entry, "system": here} if entry["system"] == "1L2Y" else entry
        for entry in config["systems"]]}
    path = _write(tmp_path / "s.yml", {**config, "output": str(tmp_path / "out"),
                                       "include_phase": ["setup"]})
    assert main(["explore", "--config", str(path), "--dry-run"]) == 2
    assert says in capsys.readouterr().err
    assert main(["explore", "--config", str(path)]) == 2
    assert says in capsys.readouterr().err
    assert not (tmp_path / "out").exists()


@pytest.mark.parametrize("config, says, code", CASES, ids=IDS)
def test_check_it_refuses_it_with_the_same_code(tmp_path, config, says, code):
    from fastmdxplora.gui.config_builder import check_config_file

    path = _write(tmp_path / "s.yml", config)
    checked = check_config_file(str(path), to_run=True, base=tmp_path)
    assert checked["ok"] is False
    assert says in checked["error"]
    assert checked["refusal"]["code"] == code
    # Opened in the form, not to run, it is not judged as a run.
    assert check_config_file(str(path))["ok"] is True


def test_check_it_through_its_route(tmp_path, monkeypatch):
    """The page's own route, for what needs no folder to read it from."""
    import urllib.request

    from fastmdxplora.gui import preview
    from fastmdxplora.gui.server import start_dashboard_session

    monkeypatch.setattr(preview, "_platform", lambda requested, precision: ("CPU", precision))
    studies = tmp_path / "lab"
    studies.mkdir()
    path = _write(studies / "s.yml", {"systems": [{"id": "a", "system": "1L2YA"}]})
    session = start_dashboard_session(output=str(studies), host="127.0.0.1", port=0,
                                      home_mode=True)
    try:
        request = urllib.request.Request(
            f"http://127.0.0.1:{session.port}/api/check-config",
            data=json.dumps({"path": str(path)}).encode(),
            headers={"Content-Type": "application/json"}, method="POST")
        with urllib.request.urlopen(request) as response:
            checked = json.load(response)
    finally:
        session.stop()
    assert checked["ok"] is False and "Could not classify" in checked["error"]


@pytest.fixture
def wire(tmp_path):
    from fastmdxplora.mcp import App, Workspace
    from tests._mcp_wire import Wire

    root = tmp_path / "work"
    root.mkdir()
    wire = Wire(App(Workspace.at(root)).server())
    yield wire, root
    wire.close()


@pytest.mark.parametrize("config, says, code", CASES, ids=IDS)
def test_an_ai_app_s_check_refuses_it_with_the_same_code(wire, config, says, code):
    wire, root = wire
    result = wire.request("tools/call", {"name": "check_study", "arguments": {
        "config": yaml.safe_dump(config)}})["result"]
    said = result["content"][0]["text"]
    assert result["isError"], said
    assert said.startswith(f"Refused ({code}):") and says in said


def test_an_ai_app_s_check_reads_a_file_from_the_workspace(wire):
    wire, root = wire
    _make_traj_files(root / "made")
    text = (root / "made" / "simulation" / "topology.pdb").read_text()
    (root / "here.pdb").write_text(text, encoding="utf-8")
    result = wire.request("tools/call", {"name": "check_study", "arguments": {
        "config": "systems:\n  - system: here.pdb\n"}})["result"]
    assert not result["isError"], result["content"][0]["text"]
    # One named from the workspace that is not there is refused, saying so.
    result = wire.request("tools/call", {"name": "check_study", "arguments": {
        "config": "systems:\n  - system: gone.pdb\n"}})["result"]
    said = result["content"][0]["text"]
    assert result["isError"] and "No structure at 'gone.pdb'" in said
    assert str(root) in said


# --------------------------------------------- what is not refused, and why

def _prepared(tmp_path: Path) -> Path:
    earlier = tmp_path / "earlier"
    (earlier / "setup").mkdir(parents=True)
    for name in ("system.pdb", "system.xml", "state.xml"):
        (earlier / "setup" / name).write_text("", encoding="utf-8")
    return earlier


@pytest.mark.parametrize("block", [
    {"simulation": {"setup_from": "earlier"}},
    {"simulation": {"prepared_from": "earlier"}},
    {"include_phase": ["simulation", "analysis"]},
    {"exclude_phase": ["setup"]},
], ids=["setup_from", "prepared_from", "include_phase", "exclude_phase"])
def test_a_system_setup_will_not_read_is_not_checked(tmp_path, monkeypatch, block):
    from fastmdxplora.config.before_run import refuse_what_can_be_seen

    monkeypatch.chdir(tmp_path)
    _prepared(tmp_path)
    config = {"systems": [{"id": "a", "system": "gone.pdb"}], **block}
    refuse_what_can_be_seen(config, base=tmp_path)
    # The same, given in the system's own block.
    if "simulation" in block:
        entry = {"id": "a", "system": "gone.pdb", **block}
        refuse_what_can_be_seen({"systems": [entry]}, base=tmp_path)
    # Without it, setup would read the file, and it is not there.
    with pytest.raises(Exception, match="No structure at 'gone.pdb'"):
        refuse_what_can_be_seen({"systems": [{"id": "a", "system": "gone.pdb"}]},
                                base=tmp_path)


def test_a_dry_run_of_a_study_from_a_prepared_system_plans_it(tmp_path, monkeypatch, capsys):
    monkeypatch.chdir(tmp_path)
    _prepared(tmp_path)
    path = _write(tmp_path / "s.yml", {"systems": [{"id": "a", "system": "gone.pdb"}],
                                       "simulation": {"setup_from": "earlier"}})
    assert main(["explore", "--config", str(path), "--dry-run"]) == 0
    assert "No structure" not in capsys.readouterr().err


def test_an_analysis_s_own_default_is_not_read_as_a_selection():
    from fastmdxplora.config.before_run import refuse_what_can_be_seen, selections_in

    config = {"systems": [{"id": "a", "system": "1L2Y"}],
              "analysis": {"options": {"water_sites": {"site_selection": "auto"}}}}
    refuse_what_can_be_seen(config)
    assert selections_in(config) == []
    # A selection given there is one, and is parsed.
    config["analysis"]["options"]["water_sites"]["site_selection"] = BAD
    with pytest.raises(Exception, match="is not an MDTraj selection"):
        refuse_what_can_be_seen(config)
    # The selection options of an analysis that takes two.
    config = {"systems": [{"id": "a", "system": "1L2Y"}],
              "analysis": {"options": {"coordination_number": {"selection_a": BAD}}}}
    with pytest.raises(Exception, match="selection_a"):
        refuse_what_can_be_seen(config)


# ------------------------------------------- the second review's leftovers

@pytest.mark.parametrize("config, says, code", CASES, ids=IDS)
def test_an_ai_app_s_save_refuses_it_and_saves_nothing(wire, config, says, code):
    """save_study saved what check_study refused, and gave it a plan id."""
    wire, root = wire
    result = wire.request("tools/call", {"name": "save_study", "arguments": {
        "config": yaml.safe_dump(config), "name": "bad"}})["result"]
    said = result["content"][0]["text"]
    assert result["isError"], said
    assert says in said and code in said
    assert not (root / "bad.yml").exists()


def test_python_is_given_the_same_check():
    """`validate_config` judges the settings; the check a run makes before
    it starts is named beside it for a script to call."""
    from fastmdxplora.config import ConfigError, refuse_what_can_be_seen

    with pytest.raises(ConfigError) as caught:
        refuse_what_can_be_seen({"systems": [{"system": "1L2YA"}]})
    assert caught.value.code == "setup.input.unrecognised"


def test_check_it_judges_the_study_your_defaults_make(tmp_path):
    """Check it judged the Config without your defaults, where the run and
    check_study judge it with them: a selection your file gives that does
    not parse was "Runs." here."""
    from fastmdxplora.config.defaults_file import DEFAULTS_FILE, read_defaults
    from fastmdxplora.gui.config_builder import check_config_file

    yours = tmp_path / DEFAULTS_FILE
    yours.write_text(f"analysis:\n  select_atoms: '{BAD}'\n", encoding="utf-8")
    path = _write(tmp_path / "s.yml", {"systems": [{"id": "a", "system": "1L2Y"}]})
    checked = check_config_file(str(path), to_run=True, base=tmp_path,
                                defaults=read_defaults(yours))
    assert checked["ok"] is False, checked
    assert checked["refusal"]["code"] == SELECTION[1]
