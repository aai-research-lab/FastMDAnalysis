"""Agent Settings says what each provider sign-in can do, and claims nothing it did not.

The settings endpoints run the provider actions in the dashboard's own
process and return only status and safe explanations. These tests cover the
actions around the main flows in ``test_agent_provider_authentication.py``:
a dashboard that is not local, a provider or method that does not exist, an
API key supplied by the environment, starting an Anthropic Console sign-in,
switching ChatGPT accounts, removing a saved key while the environment still
supplies one, testing a route that was never saved, and a key the vault
refuses. Provider modules are replaced at their public functions, so nothing
leaves the machine and no terminal or browser opens.
"""

from __future__ import annotations

import pytest


@pytest.fixture
def settings(monkeypatch, tmp_path):
    monkeypatch.setenv("FASTMDXPLORA_CONFIG_DIR", str(tmp_path))
    for name in ("ANTHROPIC_API_KEY", "OPENAI_API_KEY", "FASTMDX_MODEL_API_KEY"):
        monkeypatch.delenv(name, raising=False)
    from fastmdxplora.gui import agent_panel

    return agent_panel


def _raise(error):
    def raising(*_args, **_kwargs):
        raise error
    return raising


def test_a_hosted_dashboard_offers_no_credential_actions(monkeypatch, settings):
    from fastmdxplora.agent import openai_chatgpt as oauth

    monkeypatch.setattr(oauth, "status", _raise(AssertionError("vault must not be read")))
    for action in ("status", "connect", "activate", "disconnect", "test"):
        answer = settings.auth_endpoint({"action": action, "provider": "openai",
                                         "auth_method": "chatgpt_plan"},
                                        local_auth_allowed=False)
        assert answer["ok"] is False
        assert "available only in the local dashboard" in answer["error"]


@pytest.mark.parametrize("payload,error", [
    ({"action": "status"}, "Choose a provider in Agent Settings first."),
    ({"action": "status", "provider": "nobody"}, "Choose a provider in Agent Settings first."),
    ({"action": "connect", "provider": "anthropic", "auth_method": "chatgpt_plan"},
     "That sign-in method is not supported for this provider."),
    ({"action": "connect", "provider": "compatible", "auth_method": "anthropic_console"},
     "That sign-in method is not supported for this provider."),
    ({"action": "rotate", "provider": "openai"}, "Unknown authentication action."),
])
def test_an_unknown_provider_method_or_action_is_refused(settings, payload, error):
    assert settings.auth_endpoint(payload) == {"ok": False, "error": error}


def test_a_key_from_the_environment_is_reported_as_untested(monkeypatch, settings):
    monkeypatch.setenv("OPENAI_API_KEY", "sk-synthetic-environment")
    answer = settings.auth_endpoint({"action": "status", "provider": "openai"})
    assert answer == {"ok": True, "auth_method": "api_key", "status": {
        "available": True, "connected": None,
        "message": "API key supplied by OPENAI_API_KEY; test it to verify access."}}


def test_interactive_sign_in_is_only_for_the_methods_that_have_one(settings):
    for payload in ({"provider": "openai", "auth_method": "api_key"},
                    {"provider": "anthropic", "auth_method": "api_key"}):
        answer = settings.auth_endpoint({"action": "connect", **payload})
        assert answer == {"ok": False,
                          "error": "Interactive sign-in is not supported for this selection."}


@pytest.mark.parametrize("result,expected", [
    ({"ok": True, "message": "Finish or cancel sign-in in the opened terminal."},
     {"available": True, "connected": None, "pending": True,
      "message": "Finish or cancel sign-in in the opened terminal."}),
    ({"ok": False, "message": ""},
     {"available": False, "connected": False, "pending": False,
      "message": "Anthropic Console sign-in could not start."}),
])
def test_starting_a_console_sign_in_reports_it_pending_or_not_started(
        monkeypatch, settings, result, expected):
    from fastmdxplora.agent import anthropic_cli

    monkeypatch.setattr(anthropic_cli, "connect", lambda: dict(result))
    answer = settings.auth_endpoint({"action": "connect", "provider": "anthropic",
                                     "auth_method": "anthropic_console"})
    assert answer == {"ok": result["ok"], "auth_method": "anthropic_console",
                      "status": expected}


def test_account_switching_is_only_for_chatgpt_sign_in(settings):
    answer = settings.auth_endpoint({"action": "activate", "provider": "anthropic",
                                     "auth_method": "anthropic_console"})
    assert answer == {"ok": False,
                      "error": "Saved-account switching is available for ChatGPT sign-in."}


def test_a_refused_account_switch_is_not_reported_as_connected(monkeypatch, settings):
    from fastmdxplora.agent import openai_chatgpt as oauth

    chosen = []
    monkeypatch.setattr(oauth, "activate", lambda *, client_id: chosen.append(client_id) or {
        "ok": False, "message": ""})
    monkeypatch.setattr(oauth, "list_models", _raise(AssertionError("not listed")))
    answer = settings.auth_endpoint({"action": "activate", "provider": "openai",
                                     "auth_method": "chatgpt_plan",
                                     "client_id": "oaiapp_other"})
    assert chosen == ["oaiapp_other"]
    assert answer == {"ok": False, "auth_method": "chatgpt_plan", "status": {
        "available": True, "connected": False,
        "message": "That account could not be selected."}}


def test_a_switched_account_whose_models_cannot_be_listed_is_not_connected(
        monkeypatch, settings):
    from fastmdxplora.agent import openai_chatgpt as oauth

    accounts = [{"client_id": "oaiapp_other", "email": "other@example.invalid",
                 "subject_hint": "", "connected": True, "active": True}]
    monkeypatch.setattr(oauth, "activate", lambda *, client_id: {"ok": True, "message": "ok"})
    monkeypatch.setattr(oauth, "list_models", _raise(oauth.SignInError(
        "ChatGPT plan routing is temporarily unavailable. Try again shortly.")))
    monkeypatch.setattr(oauth, "status", lambda: {"accounts": accounts})
    answer = settings.auth_endpoint({"action": "activate", "provider": "openai",
                                     "auth_method": "chatgpt_plan",
                                     "client_id": "oaiapp_other"})
    assert answer == {"ok": False, "auth_method": "chatgpt_plan",
                      "status": {"available": True, "connected": False,
                                 "message": "ChatGPT plan routing is temporarily "
                                            "unavailable. Try again shortly."},
                      "accounts": accounts}


def test_a_switched_account_returns_its_own_models(monkeypatch, settings):
    from fastmdxplora.agent import openai_chatgpt as oauth

    options = [{"id": "gpt-plan", "label": "GPT plan"}]
    current = {"available": True, "connected": True, "message": "Signed in."}
    monkeypatch.setattr(oauth, "activate", lambda *, client_id: {"ok": True, "message": "ok"})
    monkeypatch.setattr(oauth, "list_models", lambda: options)
    monkeypatch.setattr(oauth, "status", lambda: dict(current))
    answer = settings.auth_endpoint({"action": "activate", "provider": "openai",
                                     "auth_method": "chatgpt_plan",
                                     "client_id": "oaiapp_other"})
    assert answer == {"ok": True, "auth_method": "chatgpt_plan", "status": current,
                      "model_options": options}


def test_removing_a_saved_key_says_the_environment_still_supplies_one(monkeypatch, settings):
    from fastmdxplora.agent.models import ModelChoice, _stored_api_key, save_choice

    choice = ModelChoice("openai", "gpt-5")
    save_choice(choice, key="sk-synthetic-saved")
    monkeypatch.setenv("OPENAI_API_KEY", "sk-synthetic-environment")
    answer = settings.auth_endpoint({"action": "disconnect", "provider": "openai"})
    assert answer["ok"] is True and answer["status"]["connected"] is None
    assert "OPENAI_API_KEY is still active" in answer["status"]["message"]
    assert "sk-synthetic" not in answer["status"]["message"]
    assert _stored_api_key(choice) == ""


def test_a_key_route_is_tested_only_once_it_is_saved(monkeypatch, settings):
    from fastmdxplora.agent import models

    monkeypatch.setattr(models, "probe_connection", _raise(AssertionError("not tested")))
    answer = settings.auth_endpoint({"action": "test", "provider": "openai"})
    assert answer == {"ok": False, "auth_method": "api_key", "status": {
        "available": True, "connected": False,
        "message": "Save this provider and authentication method before testing it."}}
    models.save_choice(models.ModelChoice("anthropic", "claude-sonnet-4-6"))
    assert settings.auth_endpoint({"action": "test", "provider": "openai"})["ok"] is False


def test_a_key_the_vault_refuses_is_not_saved_and_says_why(monkeypatch, settings, tmp_path):
    from fastmdxplora.agent import models

    monkeypatch.setattr(models, "_secure_keyring", lambda: None)
    answer = settings.model_endpoint({"provider": "openai", "model": "gpt-5",
                                      "api_key": "sk-synthetic-unsaved"})
    assert answer["ok"] is False
    assert answer["code"] == "environment.credentials.secure_store_unavailable"
    assert "sk-synthetic-unsaved" not in answer["error"]
    assert not (tmp_path / "model.json").exists()


def test_removing_a_saved_key_says_whether_there_was_one(settings):
    from fastmdxplora.agent.models import ModelChoice, save_choice

    save_choice(ModelChoice("openai", "gpt-5"), key="sk-synthetic-saved")
    removed = settings.auth_endpoint({"action": "disconnect", "provider": "openai"})
    assert removed["status"] == {"available": True, "connected": False,
                                 "message": "Saved API key removed from the OS credential vault."}
    again = settings.auth_endpoint({"action": "disconnect", "provider": "openai"})
    assert again["ok"] is True and again["status"]["message"] == (
        "No saved API key was found in the OS credential vault.")
