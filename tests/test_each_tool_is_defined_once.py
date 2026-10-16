"""Each tool offered to an AI model is defined once, and both interfaces offer it.

The Agent and `fastmdx mcp` each kept their own list of tools, and drifted:
an AI app could list a study's saved views and the machines a study can be
sent to, which the Agent could not; the Agent could read a file of a study
and the docs, which an AI app could not; and each read a study's record and
checked a config its own way, so the two could be told different things
about one study (the Agent was given no progress of a running study; an AI
app's check read your defaults from the workspace's top while a run reads
them from the study's folder up). Now both are built from one list
(`fastmdxplora.offered_tools`), and a tool one does not offer says why.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest
import yaml

from fastmdxplora.agent.conversation import REPLY_TOOLS
from fastmdxplora.agent.tools import LOOK_WORDS, Toolbox
from fastmdxplora.mcp.tools import TOOLS, Context, ToolError
from fastmdxplora.mcp.workspace import Workspace
from fastmdxplora.offered_tools import (
    AGENT,
    AI_APP,
    OFFERED,
    REPLIES,
    offered_by,
    offered_named,
    schema_for,
)

from tests.test_the_workspace_says_its_studies import _study


def test_every_tool_is_offered_by_both_or_kept_apart_with_its_reason():
    names = [tool.name for tool in OFFERED]
    assert len(names) == len(set(names))
    for tool in OFFERED:
        assert set(tool.offered) <= {AGENT, AI_APP} and tool.offered, tool.name
        if set(tool.offered) == {AGENT, AI_APP}:
            assert not tool.apart, tool.name
        else:
            assert len(tool.apart.split()) >= 5, f"{tool.name} is kept apart without a reason"


def test_the_agent_offers_what_the_list_gives_it_and_only_looks():
    box = Toolbox()
    offered = [tool.name for tool in offered_by(AGENT)]
    # current_view is made for each request from the page's hints.
    assert list(box.names) == [name for name in offered if name != "current_view"]
    for tool in offered_by(AGENT):
        assert tool.only_looks, f"the Agent's {tool.name} changes something"
        assert LOOK_WORDS[tool.name], tool.name
        spec = next((s for s in box.specs() if s.name == tool.name), None)
        if spec is not None:
            assert spec.parameters == schema_for(tool, AGENT)
            assert tool.what in spec.description


def test_an_ai_app_is_offered_what_the_list_gives_it():
    listed = [tool.listed() for tool in TOOLS]
    assert [t["name"] for t in listed] == [
        tool.name_for(AI_APP) for tool in offered_by(AI_APP)]
    for tool, served in zip(offered_by(AI_APP), TOOLS):
        assert served.description == tool.told_to(AI_APP)
        assert served.properties == schema_for(tool, AI_APP)["properties"]
        assert served.required == tool.required and served.acts == tool.acts
        assert served.annotations == tool.annotations


def test_the_agent_s_replies_are_named_with_why_an_ai_app_has_none():
    assert set(REPLIES) == set(REPLY_TOOLS)
    assert all(len(why.split()) >= 5 for why in REPLIES.values())
    assert not set(REPLIES) & {tool.name for tool in OFFERED}


def test_each_tool_s_function_is_there():
    from importlib import import_module

    for tool in OFFERED:
        for does in filter(None, (tool.does, tool.ai_app_does)):
            module, _, name = does.partition(":")
            assert callable(getattr(import_module(module), name)), does


def test_a_tool_the_agent_offers_too_is_the_agent_s_own_look(monkeypatch):
    """Not a second copy: what an AI app's call runs is the Agent's function."""
    from fastmdxplora.agent import tools

    seen = []
    name, (said, what, _) = "list_studies", tools._TOOLS["list_studies"]
    monkeypatch.setitem(tools._TOOLS, name,
                        (said, what, lambda box, asked: seen.append(box.for_ai_app) or "x"))
    served = next(t for t in TOOLS if t.name == name)
    assert served.run(Context(Workspace.at(Path(__file__).parent)), {}) == "x"
    assert seen == [True]


@pytest.fixture
def workspace(tmp_path):
    _study(tmp_path / "ubiquitin", means={"rmsd": (0.1234, 0.0056), "rg": (1.18, 0.002)})
    return tmp_path


def _ai_app(workspace, name, **asked):
    tool = next(t for t in TOOLS if t.name == name)
    return tool.run(Context(Workspace.at(workspace)), tool.checked(asked))


def _agent(workspace, name, **asked):
    look = Toolbox(workspace=Workspace.at(workspace)).use(name, asked)
    assert look.ok, look.said
    return look.said


def test_a_study_s_record_is_one_for_both(workspace):
    """Each its own way, the Agent's had no word on where a study stands,
    and an AI app's could not be asked for one analysis first."""
    ai_app = _ai_app(workspace, "read_study", study="ubiquitin", analysis="rg")
    agent = _agent(workspace, "read_study", study="ubiquitin", analysis="rg")
    assert ai_app == agent
    assert agent.startswith("The study at ubiquitin\nstatus: completed")
    found = agent.split("what the analyses found:\n", 1)[1]
    assert found.startswith("  rg: mean 1.1800 ± 0.0020 nm (s.e.)")


def test_a_config_is_checked_the_same_way_your_defaults_from_the_study_s_folder_up(workspace):
    """An AI app's check read your defaults from the workspace's top, a run
    from the study's folder up: a nearer file applied to the run and not
    to the plan checked. The Agent's check filled none."""
    project = workspace / "project"
    project.mkdir()
    (project / "fastmdx-defaults.yml").write_text(yaml.safe_dump(
        {"simulation": {"temperature_K": 310}}), encoding="utf-8")
    config = {"systems": [{"system": "1UBQ"}], "output": "project/ubq_310"}
    ai_app = _ai_app(workspace, "check_study", config=yaml.safe_dump(config))
    agent = _agent(workspace, "check_config", config=config)
    filled = "Your defaults (fastmdx-defaults.yml) fill what it leaves unset: " \
             "simulation.temperature_K."
    assert filled in ai_app and filled in agent
    assert agent.startswith("Accepted by the validator. The plan:")
    assert ai_app.startswith(agent)
    assert ai_app[len(agent):].strip().startswith("To run it, save it first with save_study.")
    elsewhere = _agent(workspace, "check_config", config={"systems": [{"system": "1UBQ"}]})
    assert "Your defaults" not in elsewhere


def test_the_agent_looks_at_what_an_ai_app_could(workspace, monkeypatch):
    from fastmdxplora.remote import api

    monkeypatch.setattr(api, "machines", lambda **_: [])
    monkeypatch.setattr(api, "jobs", lambda **_: [])
    said = _agent(workspace, "views_of_study", study="ubiquitin")
    assert said.startswith("Views saved with ubiquitin in the GUI: 0.")
    assert "write_scene" not in said
    machines = _agent(workspace, "list_machines")
    assert machines.startswith("No machine has been inspected from this computer.")
    assert "start_study" not in machines
    status = _agent(workspace, "remote_status")
    assert status.startswith("No job sent from here is waiting or running.")


def test_an_ai_app_reads_what_the_agent_could(workspace):
    said = _ai_app(workspace, "read_study_file", study="ubiquitin", file="manifest.json")
    assert said.startswith("manifest.json in the study at ubiquitin")
    assert json.loads(said.split(":\n", 1)[1])["phases"][0]["name"] == "setup"
    with pytest.raises(ToolError, match="outside the workspace"):
        _ai_app(workspace, "read_study_file", study="../elsewhere", file="manifest.json")
    pages = _ai_app(workspace, "read_docs")
    assert "agent" in pages


def test_a_tool_kept_apart_is_named_by_the_list():
    for name in ("save_study", "start_study", "write_scene", "ask_agent"):
        tool = offered_named(name)
        assert tool.offered == (AI_APP,) and tool.apart
    assert offered_named("check_study") is offered_named("check_config")
    assert offered_named("current_view").offered == (AGENT,)


def _defaults(folder, kelvin):
    folder.mkdir(parents=True, exist_ok=True)
    (folder / "fastmdx-defaults.yml").write_text(yaml.safe_dump(
        {"simulation": {"temperature_K": kelvin}}), encoding="utf-8")


def test_your_defaults_outside_the_workspace_are_neither_read_nor_said(workspace, tmp_path):
    """Looked for beside a study's output, a file of defaults outside the
    workspace was read, and its folder and values said, before the output
    was refused."""
    outside = tmp_path.parent / f"{tmp_path.name}-elsewhere"
    _defaults(outside / "lab", 311.5)
    config = {"systems": [{"system": "1UBQ"}], "output": str(outside / "lab" / "run")}
    with pytest.raises(ToolError) as refused:
        _ai_app(workspace, "check_study", config=yaml.safe_dump(config))
    assert "outside the workspace" in str(refused.value)
    assert "311" not in str(refused.value) and "fastmdx-defaults" not in str(refused.value)
    held = Workspace.at(workspace)
    look = Toolbox(path_for=held.path_for, workspace=held).use(
        "check_config", {"config": config})
    assert not look.ok and "outside the workspace" in look.said and "311" not in look.said


def test_a_hosted_gui_s_agent_reaches_no_other_machine(workspace, monkeypatch):
    from fastmdxplora.remote import api

    held = Workspace.at(workspace)
    hosted = Toolbox(path_for=held.path_for, workspace=held, hosted=True)
    for name in ("list_machines", "remote_status"):
        assert name not in hosted.names and f"`{name}`" not in hosted.describe()
        look = hosted.use(name, {})
        assert not look.ok and "There is no tool called" in look.said, name
    # The Agent an AI app asks (ask_agent) is held to the workspace, and is
    # not a hosted GUI: it looks at the machines.
    monkeypatch.setattr(api, "machines", lambda **_: [])
    monkeypatch.setattr(api, "jobs", lambda **_: [])
    asked = Toolbox(path_for=held.path_for, workspace=held).use("list_machines", {})
    assert asked.ok and asked.said.startswith("No machine has been inspected")


def test_a_study_sent_elsewhere_is_filled_as_it_was_checked(workspace):
    """A send read your defaults beside its config file alone; the plan
    checked read them from the study's folder up."""
    from fastmdxplora.remote.send import _with_your_defaults

    _defaults(workspace, 300)
    _defaults(workspace / "project", 310)
    config = {"systems": [{"system": "1UBQ"}], "output": "project/ubq_310"}
    (workspace / "ubq.yml").write_text(yaml.safe_dump(config), encoding="utf-8")
    sent, said = _with_your_defaults(config, workspace / "ubq.yml")
    assert sent["simulation"]["temperature_K"] == 310
    checked = _ai_app(workspace, "check_study", config="ubq.yml")
    assert "  Conditions: 310 K, 1 bar\n" in checked
    # Kept in a folder of its own, with its results folder given, as an
    # AI app's start_study sends it.
    (workspace / "sub").mkdir()
    (workspace / "sub" / "ubq.yml").write_text(yaml.safe_dump(config), encoding="utf-8")
    sent, _ = _with_your_defaults(config, workspace / "sub" / "ubq.yml",
                                  str(workspace / "project" / "ubq_310"))
    assert sent["simulation"]["temperature_K"] == 310
    # Its results folder the workspace's `run`, not one beside it: neither
    # the check nor the send takes the defaults kept beside the file.
    _defaults(workspace / "sub", 320)
    beside = {"systems": [{"system": "1UBQ"}], "output": "run"}
    (workspace / "sub" / "a.yml").write_text(yaml.safe_dump(beside), encoding="utf-8")
    checked = _ai_app(workspace, "check_study", config="sub/a.yml")
    assert "  Conditions: 300 K, 1 bar\n" in checked
    sent, _ = _with_your_defaults(beside, workspace / "sub" / "a.yml", str(workspace / "run"))
    assert sent["simulation"]["temperature_K"] == 300


def test_a_study_that_has_run_is_not_filled_and_is_said_so(workspace):
    _defaults(workspace, 300)
    config = {"systems": [{"system": "1UBQ"}], "output": "ubiquitin"}
    said = _agent(workspace, "check_config", config=config)
    assert ("Your defaults are not filled in: ubiquitin has run, and its record holds "
            "the values it ran with.") in said
    assert "fill what it leaves unset" not in said


#: Tools an AI app is offered that only look, and the Agent is not: none.
#: A look kept from the Agent is listed here with why, not only in the list.
AI_APP_ONLY_LOOKS: dict[str, str] = {}


def test_a_look_an_ai_app_is_offered_is_the_agent_s_too():
    for tool in offered_by(AI_APP):
        if tool.only_looks and AGENT not in tool.offered:
            assert tool.name in AI_APP_ONLY_LOOKS, (
                f"{tool.name} only looks: offer it to the Agent too, or say here why not")


def test_the_one_record_says_the_agent_s_words_as_they_apply_to_an_ai_app(
        workspace, monkeypatch):
    """The same record, with a fix to run: the Agent is told what its own
    `DO: run the fix` runs, an AI app what can be run here."""
    from types import SimpleNamespace

    from fastmdxplora.workspace_studies import _SAID_TO_AN_AI_APP

    fix = SimpleNamespace(argv=("resume", "ubiquitin"), decision=False,
                          as_text=lambda: "fastmdx resume ubiquitin; about 15 min here")
    monkeypatch.setattr("fastmdxplora.remedies.remedies_of", lambda here: [fix])
    agent = _agent(workspace, "read_study", study="ubiquitin")
    ai_app = _ai_app(workspace, "read_study", study="ubiquitin")
    theirs, ours = _SAID_TO_AN_AI_APP[0]
    assert theirs in agent and theirs not in ai_app and ours in ai_app
    assert ai_app.replace(ours, theirs) == agent
