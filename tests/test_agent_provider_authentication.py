"""Provider authentication remains isolated from UI and ordinary settings."""

from __future__ import annotations

import json
import urllib.error
import urllib.parse
import urllib.request
from io import BytesIO
from pathlib import Path
import hashlib
import time

import pytest


@pytest.fixture
def chatgpt_vault(monkeypatch, tmp_path):
    from fastmdxplora.agent import models
    from fastmdxplora.agent import openai_chatgpt as oauth

    values = {}

    class Vault:
        def get_password(self, service, username):
            return values.get((service, username))

        def set_password(self, service, username, value):
            values[(service, username)] = value

        def delete_password(self, service, username):
            if (service, username) not in values:
                raise ValueError("Password not found")
            del values[(service, username)]

    vault = Vault()
    monkeypatch.setattr(models, "_secure_keyring", lambda: vault)
    monkeypatch.setattr(oauth, "user_config_dir", lambda: tmp_path)
    monkeypatch.setattr(oauth, "_pending", {})
    monkeypatch.setattr(oauth, "_pending_note", "")
    return oauth, vault


def _account(*, expires_at, scopes=None):
    return {
        "client_id": "oaiapp_test",
        "subject": "subject-test",
        "issuer": "https://auth.openai.com",
        "email": "test@example.invalid",
        "ext_agent_host_id": "urn:uuid:test-host",
        "id_token": "identity-token-test",
        "access_token": "access-token-test",
        "refresh_token": "refresh-token-test",
        "token_type": "Bearer",
        "expires_at": expires_at,
        "scopes": scopes or ["openid", "offline_access", "resource.invoke",
                              "chatgpt.tokens.use.direct"],
    }


def _enable_account(oauth, account):
    oauth._write_store({"host_id": "urn:uuid:test-host",
                        "active_client_id": account["client_id"],
                        "accounts": [account]})


def test_settings_save_auth_mode_without_credentials_and_use_matching_endpoint(
        monkeypatch, tmp_path):
    from fastmdxplora.agent.models import ModelChoice, load_choice, save_choice
    from fastmdxplora.gui.agent_panel import model_endpoint

    monkeypatch.setenv("FASTMDXPLORA_CONFIG_DIR", str(tmp_path))
    choice = ModelChoice("openai", "gpt-6", auth_method="chatgpt_plan")
    assert choice.url == "https://api.openai.com/v1/responses"
    save_choice(choice)
    saved = load_choice()
    assert (saved.provider, saved.model, saved.auth_method) == (
        "openai", "gpt-6", "chatgpt_plan")
    record = json.loads(Path(tmp_path / "model.json").read_text())
    assert record == {"provider": "openai", "model": "gpt-6",
                      "auth_method": "chatgpt_plan"}

    providers = model_endpoint({})["providers"]
    openai = next(row for row in providers if row["id"] == "openai")
    siwc = next(row for row in openai["auth_methods"]
                if row["id"] == "chatgpt_plan")
    assert siwc["available"] is True
    assert siwc["requires_api_key"] is False


def test_console_login_retry_uses_pending_registration_without_assuming_saved_account(
        monkeypatch, chatgpt_vault):
    oauth, _vault = chatgpt_vault
    oauth._write_store({"host_id": "urn:uuid:retry-host", "active_client_id": "",
                        "pending_client_id": "oaiapp_pending", "accounts": []})
    monkeypatch.setattr(oauth, "_exchange", lambda _form: {
        "access_token": "access-token-retry", "refresh_token": "refresh-token-retry",
        "id_token": "signed-identity-token", "token_type": "Bearer", "expires_in": 3600,
        "scope": "openid profile email offline_access resource.invoke "
                 "chatgpt.tokens.use.direct"})
    monkeypatch.setattr(oauth, "_verify_id_token", lambda *_a, **_kw: {
        "iss": oauth.ISSUER, "sub": "retry-account", "email": "retry@example.invalid"})

    def open_authorize(url, **_kwargs):
        query = urllib.parse.parse_qs(urllib.parse.urlsplit(url).query)
        assert query["client_id"] == ["oaiapp_pending"]
        assert query["prompt"] == ["consent"]
        assert "agent_name_hint" not in query
        assert "id_token_hint" not in query
        assert "login_hint" not in query
        callback = urllib.parse.urlsplit(query["redirect_uri"][0])
        callback_query = urllib.parse.urlencode({"code": "retry-code",
                                                 "state": query["state"][0],
                                                 "client_id": "oaiapp_pending"})
        with urllib.request.urlopen(urllib.parse.urlunsplit(
                (callback.scheme, callback.netloc, callback.path, callback_query, ""))) as response:
            assert response.status == 200
        return True

    monkeypatch.setattr(oauth.webbrowser, "open", open_authorize)
    assert oauth.connect()["ok"] is True
    saved = oauth._read_store()
    assert saved["pending_client_id"] == ""
    assert saved["active_client_id"] == "oaiapp_pending"
    assert saved["accounts"][0]["subject"] == "retry-account"


def test_chatgpt_registration_callback_uses_pkce_and_keeps_tokens_out_of_settings(
        monkeypatch, tmp_path, chatgpt_vault):
    oauth, _vault = chatgpt_vault
    exchanged = []
    verified = []

    def exchange(form):
        exchanged.append(form)
        return {"access_token": "access-token-test", "refresh_token": "refresh-token-test",
                "id_token": "signed-identity-token", "token_type": "Bearer",
                "expires_in": 3600,
                "scope": "openid profile email offline_access resource.invoke "
                         "chatgpt.tokens.use.direct"}

    def verify(token, *, client_id, nonce, timeout=15.0):
        verified.append((token, client_id, nonce))
        return {"iss": oauth.ISSUER, "sub": "account-subject", "email": "a@example.invalid"}

    monkeypatch.setattr(oauth, "_exchange", exchange)
    monkeypatch.setattr(oauth, "_verify_id_token", verify)
    seen = {}

    def open_authorize(url, **_kwargs):
        parsed = urllib.parse.urlsplit(url)
        query = urllib.parse.parse_qs(parsed.query)
        seen.update(query)
        assert query["client_id"] == ["dynamic_agent_client"]
        assert query["agent_name_hint"] == ["FastMDXplora"]
        assert query["redirect_uri"][0].startswith("http://127.0.0.1:")
        assert urllib.parse.urlsplit(query["redirect_uri"][0]).path == "/auth/callback"
        assert query["resource"] == [oauth.RESOURCE]
        assert query["code_challenge_method"] == ["S256"]
        assert "chatgpt.tokens.use.direct" in query["scope"][0]
        callback = urllib.parse.urlsplit(query["redirect_uri"][0])
        callback_query = urllib.parse.urlencode({"code": "one-time-code",
                                                 "state": query["state"][0],
                                                 "client_id": "oaiapp_test"})
        with urllib.request.urlopen(urllib.parse.urlunsplit(
                (callback.scheme, callback.netloc, callback.path, callback_query, ""))) as response:
            assert response.status == 200
        return True

    monkeypatch.setattr(oauth.webbrowser, "open", open_authorize)
    result = oauth.connect()
    assert result["ok"] is True
    assert "code_verifier" in exchanged[0]
    assert exchanged[0]["client_id"] == "oaiapp_test"
    assert exchanged[0]["resource"] == oauth.RESOURCE
    assert "client_secret" not in exchanged[0]
    assert verified[0][0] == "signed-identity-token"
    assert verified[0][1] == "oaiapp_test"
    record = oauth._read_store()
    account = record["accounts"][0]
    assert record["active_client_id"] == "oaiapp_test"
    assert account["subject"] == "account-subject"
    assert oauth.access_token() == "access-token-test"
    settings = tmp_path / "model.json"
    settings.write_text(json.dumps({"provider": "openai", "model": "gpt-6",
                                    "auth_method": "chatgpt_plan"}))
    assert "access-token-test" not in settings.read_text()
    assert "refresh-token-test" not in settings.read_text()


def test_chatgpt_id_token_signature_issuer_audience_expiry_and_nonce_are_checked(
        monkeypatch):
    from fastmdxplora.agent import openai_chatgpt as oauth

    modulus = int(
        "9d0b9558969c6516a0614e724784cd45aab568f23cfc7f408f0f02fbbefee9da"
        "9d6ba0aa067f6d1ab0e28699f77810791b4840edb80614873b532c359f98b480"
        "6b9aac9b12a94150300b393fd4a872f6703d19e0fd09faa7fa1e528b52b241d"
        "ab4171ed17219a208fa0178652917c583708d344e248eb1c49233016c3ab66503", 16)
    private = int(
        "419616d71a7f525072c670a6956bc6f6b36673dfc6cc84b0878e736683cb9e0a"
        "a8f7cbef762b3d7b8a27f0a2d891120021feb8ef0e8dab8116ab93a9097012e"
        "16db7393e6f0acbb6479fded6ac595ae6a9704da188f9d8cafab4a5cd1c1620"
        "c65a968435aeeb1fceac4f249ef26aa7d59afeb999fff2b3102b0bbb9cb4f40941", 16)
    exponent = 65537
    nonce = "nonce-from-the-authorize-request"
    key_id = "test-rsa-signing-key"
    jwk = {"kid": key_id, "kty": "RSA", "use": "sig", "alg": "RS256",
           "n": oauth._b64url(modulus.to_bytes((modulus.bit_length() + 7) // 8, "big")),
           "e": oauth._b64url(exponent.to_bytes(3, "big"))}
    monkeypatch.setattr(oauth, "_fetch_json", lambda url, **_kwargs: (
        {"issuer": oauth.ISSUER, "jwks_uri": oauth.ISSUER + "/.well-known/jwks.json"}
        if url == oauth.DISCOVERY_URL else {"keys": [jwk]}))

    def signed_token(*, token_nonce=nonce, audience="oaiapp_test", expiry=None):
        header = oauth._b64url(json.dumps({"alg": "RS256", "kid": key_id},
                                          separators=(",", ":")).encode())
        claims = oauth._b64url(json.dumps({
            "iss": oauth.ISSUER, "sub": "verified-subject", "aud": audience,
            "exp": expiry if expiry is not None else int(time.time()) + 3600,
            "nonce": token_nonce,
        }, separators=(",", ":")).encode())
        encoded = f"{header}.{claims}"
        digest = hashlib.sha256(encoded.encode("ascii")).digest()
        digest_info = bytes.fromhex("3031300d060960864801650304020105000420") + digest
        size = (modulus.bit_length() + 7) // 8
        padding = size - len(digest_info) - 3
        block = b"\x00\x01" + b"\xff" * padding + b"\x00" + digest_info
        signature = pow(int.from_bytes(block, "big"), private, modulus).to_bytes(size, "big")
        return encoded + "." + oauth._b64url(signature)

    token = signed_token()
    claims = oauth._verify_id_token(token, client_id="oaiapp_test", nonce=nonce)
    assert claims["sub"] == "verified-subject"
    for bad_token, client_id, expected_nonce in (
            (signed_token(token_nonce="wrong"), "oaiapp_test", nonce),
            (signed_token(audience="different-client"), "oaiapp_test", nonce),
            (signed_token(expiry=int(time.time()) - 1), "oaiapp_test", nonce),
            (token.rsplit(".", 1)[0] + ".invalid-signature", "oaiapp_test", nonce)):
        with pytest.raises(oauth.SignInError):
            oauth._verify_id_token(bad_token, client_id=client_id, nonce=expected_nonce)


def test_chatgpt_cancelled_browser_consent_does_not_exchange_or_replace_account(
        monkeypatch, chatgpt_vault):
    oauth, _vault = chatgpt_vault
    called = []

    def open_authorize(url, **_kwargs):
        query = urllib.parse.parse_qs(urllib.parse.urlsplit(url).query)
        callback = urllib.parse.urlsplit(query["redirect_uri"][0])
        callback_query = urllib.parse.urlencode({"error": "access_denied",
                                                 "state": query["state"][0]})
        with pytest.raises(urllib.error.HTTPError) as caught:
            urllib.request.urlopen(urllib.parse.urlunsplit(
                (callback.scheme, callback.netloc, callback.path, callback_query, "")))
        assert caught.value.code == 400
        return True

    monkeypatch.setattr(oauth, "_exchange", lambda _form: called.append(True))
    monkeypatch.setattr(oauth.webbrowser, "open", open_authorize)
    result = oauth.connect()
    assert result["ok"] is True
    assert called == []
    assert not oauth._read_store()["accounts"]
    status = oauth.status()
    assert status["connected"] is False
    assert "cancelled" in status["message"].lower()


def test_invalid_identity_callback_keeps_only_the_unverified_registration_for_retry(
        monkeypatch, chatgpt_vault):
    oauth, _vault = chatgpt_vault

    monkeypatch.setattr(oauth, "_exchange", lambda _form: {
        "access_token": "access-token-test", "refresh_token": "refresh-token-test",
        "id_token": "bad-id-token", "token_type": "Bearer", "expires_in": 3600,
        "scope": "openid offline_access resource.invoke chatgpt.tokens.use.direct"})
    monkeypatch.setattr(oauth, "_verify_id_token",
                        lambda *_a, **_kw: (_ for _ in ()).throw(
                            oauth.SignInError("OpenAI identity did not validate.")))

    def open_authorize(url, **_kwargs):
        query = urllib.parse.parse_qs(urllib.parse.urlsplit(url).query)
        callback = urllib.parse.urlsplit(query["redirect_uri"][0])
        callback_query = urllib.parse.urlencode({"code": "code", "state": query["state"][0],
                                                 "client_id": "oaiapp_pending"})
        with pytest.raises(urllib.error.HTTPError) as caught:
            urllib.request.urlopen(urllib.parse.urlunsplit(
                (callback.scheme, callback.netloc, callback.path, callback_query, "")))
        assert caught.value.code == 400
        return True

    monkeypatch.setattr(oauth.webbrowser, "open", open_authorize)
    oauth.connect()
    stored = oauth._read_store()
    assert stored["pending_client_id"] == "oaiapp_pending"
    assert stored["active_client_id"] == ""
    assert stored["accounts"] == []
    assert oauth.status()["connected"] is False


def test_expired_chatgpt_access_token_refresh_rotates_tokens_under_the_saved_client(
        monkeypatch, chatgpt_vault):
    oauth, _vault = chatgpt_vault
    _enable_account(oauth, _account(expires_at=0))
    asked = []

    def exchange(form):
        asked.append(form)
        return {"access_token": "access-token-rotated",
                "refresh_token": "refresh-token-rotated", "token_type": "Bearer",
                "expires_in": 3600,
                "scope": "openid offline_access resource.invoke chatgpt.tokens.use.direct"}

    monkeypatch.setattr(oauth, "_exchange", exchange)
    assert oauth.access_token() == "access-token-rotated"
    assert asked == [{"grant_type": "refresh_token", "client_id": "oaiapp_test",
                      "refresh_token": "refresh-token-test", "resource": oauth.RESOURCE}]
    saved = oauth._read_store()["accounts"][0]
    assert saved["refresh_token"] == "refresh-token-rotated"
    assert saved["access_token"] == "access-token-rotated"


def test_expired_revoked_refresh_grant_is_cleared_but_registration_remains(
        monkeypatch, chatgpt_vault):
    oauth, _vault = chatgpt_vault
    _enable_account(oauth, _account(expires_at=0))
    monkeypatch.setattr(oauth, "_exchange", lambda _form: (_ for _ in ()).throw(
        oauth.InvalidGrant("The refresh token expired.")))
    with pytest.raises(oauth.InvalidGrant):
        oauth.access_token()
    saved = oauth._read_store()
    assert saved["active_client_id"] == ""
    assert saved["accounts"][0]["client_id"] == "oaiapp_test"
    assert not any(name in saved["accounts"][0] for name in
                   ("access_token", "refresh_token", "id_token", "scopes"))


def test_chatgpt_disconnect_revokes_and_clears_tokens_but_keeps_client_registration(
        monkeypatch, chatgpt_vault):
    oauth, _vault = chatgpt_vault
    _enable_account(oauth, _account(expires_at=9999999999))
    monkeypatch.setattr(oauth, "_fetch_json", lambda _url: {
        "revocation_endpoint": oauth.ISSUER + "/api/accounts/oauth/revoke"})
    sent = []

    class Reply(BytesIO):
        status = 200

        def __enter__(self):
            return self

        def __exit__(self, *_args):
            self.close()

    def urlopen(request, timeout=None):
        sent.append(request)
        return Reply(b"")

    monkeypatch.setattr(urllib.request, "urlopen", urlopen)
    result = oauth.disconnect()
    assert result["remote_revoked"] is True
    body = urllib.parse.parse_qs(sent[0].data.decode("ascii"))
    assert body == {"token": ["refresh-token-test"],
                    "token_type_hint": ["refresh_token"],
                    "client_id": ["oaiapp_test"]}
    saved = oauth._read_store()
    assert saved["active_client_id"] == ""
    assert saved["accounts"][0]["client_id"] == "oaiapp_test"
    assert not any(name in saved["accounts"][0] for name in
                   ("access_token", "refresh_token", "id_token", "scopes"))


def test_chatgpt_model_catalog_uses_public_list_and_only_visible_models(
        monkeypatch, chatgpt_vault):
    oauth, _vault = chatgpt_vault
    monkeypatch.setattr(oauth, "access_token", lambda: "route-test-token")

    class Reply(BytesIO):
        def __enter__(self):
            return self

        def __exit__(self, *_args):
            self.close()

    requests = []

    def urlopen(request, timeout=None):
        requests.append(request)
        return Reply(json.dumps({"models": [
            {"slug": "gpt-visible", "display_name": "GPT Visible", "visibility": "list"},
            {"slug": "gpt-hidden", "display_name": "GPT Hidden", "visibility": "private"},
        ]}).encode())

    monkeypatch.setattr(urllib.request, "urlopen", urlopen)
    assert oauth.list_models() == [{"id": "gpt-visible", "label": "GPT Visible"}]
    assert requests[0].full_url == "https://api.openai.com/v1/models"
    assert requests[0].get_header("Authorization") == "Bearer route-test-token"


def test_responses_route_sends_full_context_and_reads_function_calls(
        monkeypatch):
    from fastmdxplora.agent import turns
    from fastmdxplora.agent.turns import ToolCall, ToolSpec

    events = [
        {"type": "response.output_text.delta", "delta": "Checking."},
        {"type": "response.output_item.added", "item": {
            "id": "item_1", "type": "function_call", "call_id": "call_1",
            "name": "read_current_view", "arguments": ""}},
        {"type": "response.function_call_arguments.delta", "item_id": "item_1",
         "delta": "{\"selected\":true}"},
        {"type": "response.output_item.done", "item": {
            "id": "item_1", "type": "function_call", "call_id": "call_1",
            "name": "read_current_view", "arguments": "{\"selected\":true}"}},
        {"type": "response.completed", "response": {
            "status": "completed", "usage": {"input_tokens": 18, "output_tokens": 4},
            "output": [{"id": "item_1", "type": "function_call",
                        "call_id": "call_1", "name": "read_current_view",
                        "arguments": "{\"selected\":true}"}]}},
    ]

    class Reply:
        def __enter__(self):
            return self

        def __exit__(self, *_args):
            return None

        def __iter__(self):
            return iter([("data: " + json.dumps(event) + "\n\n").encode()
                         for event in events])

    sent = []

    def urlopen(request, timeout=None):
        sent.append(request)
        return Reply()

    monkeypatch.setattr(urllib.request, "urlopen", urlopen)
    text = []
    called = []
    turn = turns.take_turn(
        url="https://api.openai.com/v1/responses",
        headers={"Authorization": "Bearer route-test-token"}, shape="responses",
        model="gpt-visible", system="Use verified study evidence.",
        messages=[{"role": "user", "text": "Inspect this study."},
                  {"role": "assistant", "text": "Reading context.",
                   "calls": [ToolCall("call_old", "read_current_view", {})]},
                  {"role": "results", "results": [{"id": "call_old",
                                                       "name": "read_current_view",
                                                       "content": "done"}]},
                  {"role": "user", "text": "Continue."}],
        tools=[ToolSpec("read_current_view", "Read current view", {
            "type": "object", "properties": {"selected": {"type": "boolean"}}})],
        label="OpenAI", on_text=text.append, on_call=called.append)
    body = json.loads(sent[0].data.decode("utf-8"))
    assert sent[0].full_url == "https://api.openai.com/v1/responses"
    assert body["instructions"] == "Use verified study evidence."
    assert body["input"][-2] == {"type": "function_call_output", "call_id": "call_old",
                                  "output": "done"}
    assert body["store"] is False and body["stream"] is True
    assert not ({"previous_response_id", "max_output_tokens", "tool_search"} & body.keys())
    assert body["input"][0] == {
        "type": "additional_tools", "role": "developer",
        "tools": [{"type": "function", "name": "read_current_view",
                   "description": "Read current view",
                   "parameters": {"type": "object", "properties": {
                       "selected": {"type": "boolean"}}}}],
    }
    assert turn.text == "Checking."
    assert text == ["Checking."]
    assert called == ["read_current_view"]
    assert turn.calls == (ToolCall("call_1", "read_current_view", {"selected": True}),)
    assert turn.usage.input_tokens == 18 and turn.usage.output_tokens == 4


def test_responses_stream_requires_completed_and_reports_incomplete_or_failed(monkeypatch):
    from fastmdxplora.agent import turns
    from fastmdxplora.refusals import StudyError

    class Reply:
        def __init__(self, events):
            self.events = events

        def __enter__(self):
            return self

        def __exit__(self, *_args):
            return None

        def __iter__(self):
            return iter([("data: " + json.dumps(event) + "\n\n").encode()
                         for event in self.events])

    for events in ([{"type": "response.output_text.delta", "delta": "partial"}],
                   [{"type": "response.incomplete", "response": {
                       "incomplete_details": {"reason": "max_output_tokens"}}}],
                   [{"type": "response.failed", "response": {
                       "error": {"code": "subscription_sharing_usage_limit_exceeded"}}}]):
        monkeypatch.setattr(urllib.request, "urlopen", lambda *_a, events=events, **_kw: Reply(events))
        with pytest.raises(StudyError):
            turns.take_turn(url="https://api.openai.com/v1/responses",
                            headers={"Authorization": "Bearer route-test-token"},
                            shape="responses", model="gpt-visible", system="",
                            messages=[{"role": "user", "text": "hi"}], tools=[],
                            label="OpenAI")


def test_anthropic_console_profile_never_inherits_or_reads_other_cli_credentials(
        monkeypatch, tmp_path):
    from fastmdxplora.agent import anthropic_cli

    monkeypatch.setattr(anthropic_cli, "config_dir", lambda: tmp_path / "anthropic-console")
    monkeypatch.setenv("ANTHROPIC_API_KEY", "must-not-inherit")
    monkeypatch.setenv("ANTHROPIC_CONFIG_DIR", str(tmp_path / "some-other-profile"))
    monkeypatch.setenv("ANTHROPIC_BASE_URL", "https://not-used.invalid")
    env = anthropic_cli._environment()
    assert env["ANTHROPIC_CONFIG_DIR"] == str(tmp_path / "anthropic-console")
    assert "ANTHROPIC_API_KEY" not in env
    assert "ANTHROPIC_BASE_URL" not in env
    assert str(tmp_path / "some-other-profile") not in env.values()


def test_anthropic_console_uses_documented_bearer_token_command_and_logout(monkeypatch):
    from subprocess import CompletedProcess

    from fastmdxplora.agent import anthropic_cli

    commands = []

    def run(args, **_kwargs):
        commands.append(args)
        return CompletedProcess(args, 0, stdout="console-oauth-token\n", stderr="")

    monkeypatch.setattr(anthropic_cli, "executable", lambda: "ant")
    monkeypatch.setattr(anthropic_cli, "_run", run)
    assert anthropic_cli.access_token() == "console-oauth-token"
    assert commands[-1] == ["auth", "print-credentials", "--access-token"]
    assert anthropic_cli.disconnect()["ok"] is True
    assert commands[-1] == ["auth", "logout", "--all"]


def test_unsupported_claude_subscription_sign_in_stays_unavailable_and_hosted_oauth_cannot_save(
        monkeypatch, tmp_path):
    from fastmdxplora.agent import anthropic_cli
    from fastmdxplora.agent.models import load_choice
    from fastmdxplora.gui.agent_panel import model_endpoint

    monkeypatch.setenv("FASTMDXPLORA_CONFIG_DIR", str(tmp_path))
    monkeypatch.setattr(anthropic_cli, "executable", lambda: None)
    response = model_endpoint({})
    claude = next(row for row in response["providers"] if row["id"] == "anthropic")
    console = next(row for row in claude["auth_methods"]
                   if row["id"] == "anthropic_console")
    assert console["available"] is False
    assert "Claude.ai account sign-in is not supported" in console["reason"]
    assert console["requires_api_key"] is False

    rejected = model_endpoint({"provider": "openai", "model": "gpt-5",
                               "auth_method": "chatgpt_plan"},
                              local_auth_allowed=False)
    assert rejected["ok"] is False
    assert "only in the local dashboard" in rejected["error"]
    assert load_choice() is None


def test_dashboard_test_route_reports_missing_chatgpt_sign_in_without_network_access(
        monkeypatch, tmp_path):
    from fastmdxplora.gui.agent_panel import auth_endpoint

    monkeypatch.setenv("FASTMDXPLORA_CONFIG_DIR", str(tmp_path))
    response = auth_endpoint({"action": "test", "provider": "openai",
                              "auth_method": "chatgpt_plan", "model": "gpt-5"})
    assert response["ok"] is True
    assert response["status"]["connected"] is False
    assert "No active ChatGPT sign-in" in response["status"]["message"]


def test_live_agent_settings_switch_provider_method_and_persist_the_real_selection(
        monkeypatch, tmp_path):
    pytest.importorskip("playwright.sync_api")
    from playwright.sync_api import sync_playwright

    from fastmdxplora.agent import models
    from fastmdxplora.gui.server import start_dashboard_session
    from tests.test_the_drawing_scripts_run_in_a_browser import _write_study

    monkeypatch.setattr(models, "model_path", lambda: tmp_path / "model.json")
    study = _write_study(tmp_path / "study")
    session = start_dashboard_session(output=str(study), host="127.0.0.1", port=0)
    try:
        with sync_playwright() as pw:
            browser = pw.chromium.launch(args=["--enable-unsafe-swiftshader"])
            page = browser.new_page(viewport={"width": 1440, "height": 900})
            page.set_default_timeout(30000)
            page.goto(session.url + "#agent", wait_until="domcontentloaded")
            page.wait_for_function("() => document.body.classList.contains('state-ready')")
            page.wait_for_function("() => document.getElementById('agent-provider').options.length >= 3")
            page.click("#agent-settings-open")
            page.wait_for_selector("#agent-settings:not([hidden])")

            page.locator("#agent-provider").select_option("openai")
            page.wait_for_function("() => document.getElementById('agent-model').value === 'gpt-5'")
            assert page.locator("#agent-key").is_visible()
            page.locator("#agent-auth-method").select_option("chatgpt_plan")
            assert page.locator("#agent-key").is_hidden()
            assert page.locator("#agent-auth-connect").inner_text() == "Continue with ChatGPT"
            page.wait_for_function("() => document.getElementById('agent-auth-status')"
                                   ".textContent.includes('No ChatGPT account is connected')")
            assert "No ChatGPT account is connected" in page.locator("#agent-auth-status").inner_text()

            # Exercise the dashboard's real auth test endpoint while signed out.
            page.locator("#agent-auth-test").evaluate("element => element.click()")
            page.wait_for_function("() => document.getElementById('agent-auth-status')"
                                   ".textContent.includes('No active ChatGPT sign-in')")
            page.locator("#agent-auth-disconnect").evaluate("element => element.click()")
            page.wait_for_function("() => document.getElementById('agent-auth-status')"
                                   ".textContent.includes('No saved ChatGPT sign-in')")

            page.click("#agent-save-model")
            page.wait_for_selector("#agent-settings[hidden]", state="hidden")

            page.reload(wait_until="domcontentloaded")
            page.wait_for_function("() => document.body.classList.contains('state-ready')")
            page.click("#agent-settings-open")
            page.wait_for_function("() => document.getElementById('agent-provider').value === 'openai'")
            assert page.locator("#agent-auth-method").input_value() == "chatgpt_plan"
            assert "Saved for Agent requests: OpenAI · gpt-5 · ChatGPT plan sign-in" \
                in page.locator("#agent-model-current").inner_text()

            page.locator("#agent-provider").select_option("anthropic")
            page.wait_for_function(
                "() => document.getElementById('agent-model').value === 'claude-sonnet-4-6'")
            assert page.locator("#agent-key").is_visible()
            console_option = page.locator(
                "#agent-auth-method option[value='anthropic_console']")
            from fastmdxplora.agent.anthropic_cli import can_launch_login, executable

            console_available = bool(executable() and can_launch_login())
            assert console_option.evaluate("option => option.disabled") is (not console_available)
            if console_available:
                page.locator("#agent-auth-method").select_option("anthropic_console")
                assert page.locator("#agent-key").is_hidden()
                assert page.locator("#agent-auth-connect").is_visible()
            assert "Claude.ai account sign-in is not supported" \
                in page.locator("#agent-auth-help").inner_text()
            summary = page.locator("#agent-model-current").inner_text()
            assert "Saved for Agent requests: OpenAI · gpt-5 · ChatGPT plan sign-in" in summary
            assert "Current selection (unsaved): Anthropic · claude-sonnet-4-6 · API key" in summary
            browser.close()
    finally:
        session.server.shutdown()
