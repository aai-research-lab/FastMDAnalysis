"""A signed-in account's token goes to its own provider, from this computer.

Found in review of the sign-in work:

- The settings page keeps the base URL field filled while it is hidden, and
  saved it with whatever provider was chosen. With Anthropic Console sign-in
  chosen after an OpenAI-compatible server, every Agent turn sent the
  Console token to that server, and the connection test (which asks the
  fixed endpoint) still passed.
- A GUI served beyond this computer refuses sign-in, but its settings route
  still returned the accounts signed in here (their e-mails), refreshed their
  tokens, and its Agent turns used them.
- ``ant auth print-credentials`` output was sent as a token whatever it was.
- Profiles retired by a disconnect kept the credentials a late login wrote.
- The OAuth callback's port could be shared (SO_REUSEADDR).

Every key and token here is synthetic and no request leaves the machine.
"""

from __future__ import annotations

import json
from subprocess import CompletedProcess

import pytest


@pytest.fixture
def config(monkeypatch, tmp_path):
    monkeypatch.setenv("FASTMDXPLORA_CONFIG_DIR", str(tmp_path))
    for name in ("ANTHROPIC_API_KEY", "OPENAI_API_KEY", "FASTMDX_MODEL_API_KEY"):
        monkeypatch.delenv(name, raising=False)
    return tmp_path


def _saved(config, record):
    (config / "model.json").write_text(json.dumps(record), encoding="utf-8")


@pytest.mark.parametrize("auth_method", ["anthropic_console", "api_key"])
def test_a_signed_in_token_or_a_key_is_never_sent_to_a_url_left_from_another_provider(
        config, auth_method):
    from fastmdxplora.agent.models import PROVIDERS, ModelChoice, load_choice

    # As the page saved it: the hidden URL field still held Ollama's.
    _saved(config, {"provider": "anthropic", "model": "claude-sonnet-4-6",
                    "auth_method": auth_method, "base_url": "http://192.0.2.7:11434/v1"})
    loaded = load_choice()
    assert loaded is not None and loaded.base_url == ""
    assert loaded.url == PROVIDERS["anthropic"]["url"]
    # Even given one directly, a signed-in account's token keeps its endpoint.
    given = ModelChoice("anthropic", "m", "http://192.0.2.7:11434/v1", "anthropic_console")
    assert given.url == PROVIDERS["anthropic"]["url"]
    plan = ModelChoice("openai", "m", "http://192.0.2.7:11434/v1", "chatgpt_plan")
    assert "192.0.2.7" not in plan.url


def test_saving_drops_a_base_url_for_a_provider_with_its_own_endpoint(config, monkeypatch):
    from fastmdxplora.agent import models
    from fastmdxplora.gui import agent_panel

    monkeypatch.setattr(agent_panel, "_model_options", lambda choice: [])
    monkeypatch.setattr(agent_panel, "_key_availability", lambda choice: {"available": True})
    monkeypatch.setattr(models, "_vault", lambda: None, raising=False)
    saved = agent_panel.model_endpoint({"provider": "openai", "model": "gpt-5",
                                        "auth_method": "api_key",
                                        "base_url": "http://192.0.2.7:8000/v1"})
    stored = json.loads((config / "model.json").read_text(encoding="utf-8"))
    assert saved.get("ok") is not False, saved
    assert "base_url" not in stored
    # An OpenAI-compatible server keeps its own.
    agent_panel.model_endpoint({"provider": "compatible", "model": "llama3.1",
                                "base_url": "http://localhost:11434/v1"})
    stored = json.loads((config / "model.json").read_text(encoding="utf-8"))
    assert stored["base_url"] == "http://localhost:11434/v1"


def test_a_hosted_gui_neither_shows_nor_refreshes_the_accounts_signed_in_here(
        config, monkeypatch):
    from fastmdxplora.agent import openai_chatgpt
    from fastmdxplora.gui import agent_panel

    _saved(config, {"provider": "openai", "model": "gpt-5", "auth_method": "chatgpt_plan"})
    monkeypatch.setattr(openai_chatgpt, "status",
                        lambda: pytest.fail("asked the ChatGPT accounts for a hosted page"))
    monkeypatch.setattr(agent_panel, "_model_options",
                        lambda choice: pytest.fail("listed AI models with a signed-in account"))
    read = agent_panel.model_endpoint({}, local_auth_allowed=False)
    assert read["ok"] is True
    status = read["auth_status"]
    assert status["available"] is False and "only in the local dashboard" in status["message"]
    assert "@" not in json.dumps(read)


def test_a_hosted_gui_does_not_ask_with_an_account_signed_in_here(config, monkeypatch):
    from fastmdxplora.agent import models
    from fastmdxplora.gui import agent_panel

    _saved(config, {"provider": "anthropic", "model": "claude-sonnet-4-6",
                    "auth_method": "anthropic_console"})
    monkeypatch.setattr(models, "completion_for",
                        lambda *a, **kw: pytest.fail("asked with a signed-in account"))
    import fastmdxplora.agent as agent

    monkeypatch.setattr(agent, "completion_for",
                        lambda *a, **kw: pytest.fail("asked with a signed-in account"))
    answer = agent_panel.propose_endpoint({"request": "chignolin for 2 ns"},
                                          path_for=lambda path: path)
    assert answer["ok"] is False
    assert "signed-in account" in answer["error"]
    assert answer["code"] == "environment.credentials.absent"


@pytest.mark.parametrize("printed", [
    "notice: a newer ant is available\nsk-ant-oat01-synthetic",
    "sk-ant-oat01-synthetic extra",
    "\x1b[2Ksk-ant-oat01-synthetic",
])
def test_anything_but_one_token_from_ant_is_refused_and_not_repeated(monkeypatch, printed):
    from fastmdxplora.agent import anthropic_cli

    monkeypatch.setattr(anthropic_cli, "executable", lambda: "/opt/synthetic/ant")
    monkeypatch.setattr(anthropic_cli, "_run",
                        lambda args, **kw: CompletedProcess(args, 0, stdout=printed, stderr=""))
    with pytest.raises(anthropic_cli.ConsoleSignInError) as refused:
        anthropic_cli.access_token()
    assert "synthetic" not in str(refused.value)


def test_a_relative_path_entry_does_not_find_ant(monkeypatch, tmp_path):
    import os

    from fastmdxplora.agent import anthropic_cli

    if os.name == "nt":
        pytest.skip("POSIX PATH semantics")
    planted = tmp_path / "ant"
    planted.write_text("#!/bin/sh\nexit 0\n", encoding="utf-8")
    planted.chmod(0o755)
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("PATH", "." + os.pathsep + str(tmp_path / "nowhere"))
    assert anthropic_cli.executable() is None
    monkeypatch.setenv("PATH", str(tmp_path))
    assert anthropic_cli.executable() == str(planted)


def test_profiles_retired_by_a_disconnect_are_signed_out_and_removed(config, monkeypatch):
    from fastmdxplora.agent import anthropic_cli

    current = "a" * 32
    (config / "anthropic-console-profile.json").write_text(
        json.dumps({"generation": current}), encoding="utf-8")
    keep = config / ("anthropic-console-" + current)
    retired = [config / ("anthropic-console-" + "b" * 32), config / "anthropic-console"]
    other = config / "anthropic-console-notes"
    outside = config.parent / "outside-profile"
    for folder in [keep, *retired, other, outside]:
        folder.mkdir()
        (folder / "credentials.json").write_text("{}", encoding="utf-8")
    link = config / ("anthropic-console-" + "c" * 32)
    link.symlink_to(outside, target_is_directory=True)
    signed_out = []
    monkeypatch.setattr(anthropic_cli, "executable", lambda: "/opt/synthetic/ant")
    monkeypatch.setattr(anthropic_cli, "_run", lambda args, profile=None, **kw: (
        signed_out.append((args, profile)) or CompletedProcess(args, 0, "", "")))
    anthropic_cli._clear_retired_profiles()
    assert keep.is_dir() and other.is_dir() and (outside / "credentials.json").is_file()
    assert link.is_symlink()
    assert not any(folder.exists() for folder in retired)
    assert sorted(str(profile) for _, profile in signed_out) == sorted(map(str, retired))
    assert all(args == ["auth", "logout", "--all"] for args, _ in signed_out)


def test_the_sign_in_callback_port_is_not_shared():
    import socket

    from fastmdxplora.agent.openai_chatgpt import _LoopbackServer

    assert _LoopbackServer.allow_reuse_address is False
    server = _LoopbackServer(("127.0.0.1", 0), object)
    try:
        if hasattr(socket, "SO_EXCLUSIVEADDRUSE"):
            assert server.socket.getsockopt(socket.SOL_SOCKET, socket.SO_EXCLUSIVEADDRUSE) == 1
        assert server.socket.getsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR) == 0
    finally:
        server.server_close()
