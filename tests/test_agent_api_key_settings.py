"""API-key model settings survive reloads without restoring removed auth state."""

from __future__ import annotations

import importlib
import importlib.util
import inspect
import json
import os
from urllib.error import HTTPError
from urllib.request import Request, urlopen

import pytest


class _Response:
    def __init__(self, body: bytes):
        self.body = body

    def __enter__(self):
        return self

    def __exit__(self, *_args):
        return False

    def read(self):
        return self.body


def test_an_openai_choice_and_key_are_saved_and_reloaded(tmp_path, monkeypatch):
    from fastmdxplora.agent.models import ModelChoice, _key_for, load_choice, save_choice

    monkeypatch.delenv("OPENAI_API_KEY", raising=False)
    path = tmp_path / "model.json"
    choice = ModelChoice("openai", "gpt-5")
    save_choice(choice, key="saved-test-key", path=path)

    assert load_choice(path) == choice
    assert _key_for(choice, path) == "saved-test-key"
    record = json.loads(path.read_text(encoding="utf-8"))
    assert record["provider"] == "openai"
    assert record["api_key"] == "saved-test-key"
    if os.name != "nt":
        assert path.stat().st_mode & 0o077 == 0


def test_the_environment_key_wins_over_the_saved_key(tmp_path, monkeypatch):
    from fastmdxplora.agent.models import ModelChoice, _key_for, save_choice

    path = tmp_path / "model.json"
    choice = ModelChoice("openai", "gpt-5")
    save_choice(choice, key="saved-test-key", path=path)
    monkeypatch.setenv("OPENAI_API_KEY", "environment-test-key")

    assert _key_for(choice, path) == "environment-test-key"


def test_a_saved_key_is_not_used_for_another_provider_or_server(tmp_path, monkeypatch):
    from fastmdxplora.agent.models import ModelChoice, _key_for, save_choice
    from fastmdxplora.refusals import StudyError

    monkeypatch.delenv("OPENAI_API_KEY", raising=False)
    monkeypatch.delenv("FASTMDX_MODEL_API_KEY", raising=False)
    path = tmp_path / "model.json"
    save_choice(ModelChoice("openai", "gpt-5"), key="saved-test-key", path=path)

    with pytest.raises(StudyError):
        _key_for(ModelChoice("compatible", "model-x", "https://localhost/v1"), path)
    save_choice(ModelChoice("compatible", "model-x", "https://one.example/v1"),
                key="other-test-key", path=path)
    with pytest.raises(StudyError):
        _key_for(ModelChoice("compatible", "model-x", "https://two.example/v1"), path)


def test_openai_completions_use_the_chat_completions_contract(tmp_path, monkeypatch):
    from fastmdxplora.agent.models import ModelChoice, completion_for, save_choice

    monkeypatch.delenv("OPENAI_API_KEY", raising=False)
    path = tmp_path / "model.json"
    choice = ModelChoice("openai", "gpt-5")
    save_choice(choice, key="request-test-key", path=path)
    requests = []

    def fake_urlopen(request, *, timeout):
        requests.append((request, timeout))
        return _Response(b'{"choices":[{"message":{"content":"SAY: ready"}}]}')

    monkeypatch.setattr("urllib.request.urlopen", fake_urlopen)
    assert completion_for(choice, path=path)("SAY ready") == "SAY: ready"

    request, timeout = requests[0]
    assert request.full_url == "https://api.openai.com/v1/chat/completions"
    assert request.get_method() == "POST"
    assert request.get_header("Authorization") == "Bearer request-test-key"
    assert json.loads(request.data) == {
        "model": "gpt-5",
        "messages": [{"role": "user", "content": "SAY ready"}],
    }
    assert timeout == 120.0


def test_openai_models_are_listed_from_v1_models(tmp_path, monkeypatch):
    from fastmdxplora.agent.models import ModelChoice, list_models, save_choice

    monkeypatch.delenv("OPENAI_API_KEY", raising=False)
    path = tmp_path / "model.json"
    choice = ModelChoice("openai", "gpt-5")
    save_choice(choice, key="models-test-key", path=path)
    requests = []

    def fake_urlopen(request, *, timeout):
        requests.append((request, timeout))
        return _Response(b'{"data":[{"id":"gpt-5"},{"id":"gpt-5-mini"}]}')

    monkeypatch.setattr("urllib.request.urlopen", fake_urlopen)
    assert list_models(choice, path=path) == ("gpt-5", "gpt-5-mini")

    request, timeout = requests[0]
    assert request.full_url == "https://api.openai.com/v1/models"
    assert request.get_method() == "GET"
    assert request.get_header("Authorization") == "Bearer models-test-key"
    assert timeout == 10.0


def test_the_agent_model_command_saves_an_openai_choice(tmp_path, monkeypatch, capsys):
    from fastmdxplora.agent.models import ModelChoice, _key_for, load_choice
    from fastmdxplora.cli.main import main

    monkeypatch.setenv("FASTMDXPLORA_CONFIG_DIR", str(tmp_path / "settings"))
    monkeypatch.delenv("OPENAI_API_KEY", raising=False)
    monkeypatch.setattr("fastmdxplora.agent.models.list_models", lambda *_a, **_kw: ("gpt-5",))
    answers = iter(("1", "gpt-5", "cli-test-key"))
    monkeypatch.setattr("builtins.input", lambda _prompt="": next(answers))

    assert main(["agent", "model"]) == 0
    output = capsys.readouterr().out
    choice = load_choice()
    assert choice == ModelChoice("openai", "gpt-5")
    assert _key_for(choice) == "cli-test-key"
    assert "cli-test-key" not in output


def test_stale_model_json_is_ignored_then_replaced_without_legacy_fields(
        tmp_path, monkeypatch, capsys):
    from fastmdxplora.agent.models import (
        ModelChoice,
        _key_for,
        describe_choice,
        load_choice,
    )
    from fastmdxplora.cli.main import main
    from fastmdxplora.gui.agent_panel import model_endpoint

    settings = tmp_path / "settings"
    settings.mkdir()
    model_file = settings / "model.json"
    model_file.write_text(json.dumps({
        "provider": "anthropic",
        "model": "old-model",
        "cache_for": "one-hour",
        "auth_method": "chatgpt_plan",
        "oauth_access_token": "stale-test-token",
        "account_id": "stale-account",
        "session": {"id": "stale-session"},
    }), encoding="utf-8")
    monkeypatch.setenv("FASTMDXPLORA_CONFIG_DIR", str(settings))
    monkeypatch.delenv("OPENAI_API_KEY", raising=False)
    monkeypatch.delenv("FASTMDX_MODEL_API_KEY", raising=False)

    assert load_choice() is None
    assert describe_choice() == "No AI model is set. Run `fastmdx agent model` to pick one."
    gui = model_endpoint({})
    assert gui["current"] is None
    assert gui["key_set"] is False
    assert [provider["id"] for provider in gui["providers"]] == ["openai", "compatible"]
    assert "stale-test-token" not in json.dumps(gui)

    monkeypatch.setattr("builtins.input", lambda _prompt="": (_ for _ in ()).throw(EOFError()))
    assert main(["agent", "model"]) == 1
    output = capsys.readouterr().out
    assert "No AI model is set yet." in output
    assert "Anthropic" not in output and "signed in" not in output.lower()

    # A cache field does not prevent a valid OpenAI choice from loading.
    model_file.write_text(json.dumps({
        "provider": "openai", "model": "gpt-5", "cache_for": "one-hour",
    }), encoding="utf-8")
    assert load_choice() == ModelChoice("openai", "gpt-5")

    monkeypatch.setattr("fastmdxplora.agent.models.list_model_options",
                        lambda _choice: [{"id": "gpt-5", "label": "gpt-5"}])
    saved = model_endpoint({
        "provider": "openai", "model": "gpt-5", "api_key": "new-test-key",
    })
    assert saved["ok"] is True and saved["key_set"] is True
    assert saved["current"] == {"provider": "openai", "model": "gpt-5"}
    assert "new-test-key" not in json.dumps(saved)
    assert load_choice() == ModelChoice("openai", "gpt-5")
    assert _key_for(load_choice()) == "new-test-key"
    stored = json.loads(model_file.read_text(encoding="utf-8"))
    assert set(stored) == {"provider", "model", "api_key"}
    assert "cache_for" not in stored
    assert not any(name in stored for name in (
        "auth_method", "oauth_access_token", "account_id", "session"))


def test_removed_agent_modules_and_names_cannot_be_imported():
    for name in (
        "fastmdxplora.agent.anthropic_cli",
        "fastmdxplora.agent.openai_chatgpt",
    ):
        assert importlib.util.find_spec(name) is None

    models = importlib.import_module("fastmdxplora.agent.models")
    turns = importlib.import_module("fastmdxplora.agent.turns")
    panel = importlib.import_module("fastmdxplora.gui.agent_panel")
    server = importlib.import_module("fastmdxplora.gui.server")
    assert set(models.PROVIDERS) == {"openai", "compatible"}
    for module, names in (
        (models, ("CACHE_FOR", "cache_for", "save_cache_for", "ANTHROPIC_API_KEY")),
        (turns, ("_responses_body", "_responses_reply", "CACHE_FOR")),
        (panel, ("auth_endpoint", "auth_status_endpoint", "local_auth_allowed")),
        (server, ("local_auth_allowed",)),
    ):
        for name in names:
            assert not hasattr(module, name), f"{module.__name__}.{name} remains"
    assert "cache_for" not in inspect.signature(turns.take_turn).parameters


def test_removed_auth_routes_are_unknown_and_not_public(tmp_path, monkeypatch):
    from fastmdxplora.gui.server import (
        POSTS_ANSWERED_BEYOND_LOOPBACK,
        start_dashboard_session,
    )

    monkeypatch.setenv("FASTMDXPLORA_CONFIG_DIR", str(tmp_path / "settings"))
    session = start_dashboard_session(output=str(tmp_path), host="127.0.0.1", port=0)

    def post(path: str):
        request = Request(session.url.rstrip("/") + path, data=b"{}", method="POST",
                          headers={"content-type": "application/json"})
        try:
            with urlopen(request, timeout=10) as response:
                return response.status, response.read()
        except HTTPError as exc:
            return exc.code, exc.read()

    try:
        unknown = post("/api/route-that-does-not-exist")
        for path in ("/api/agent/auth",):
            assert path not in POSTS_ANSWERED_BEYOND_LOOPBACK
            assert post(path) == unknown
    finally:
        session.server.shutdown()
