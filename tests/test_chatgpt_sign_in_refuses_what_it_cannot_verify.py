"""Sign in with ChatGPT refuses every response, record and race it cannot verify.

The happy path and the callback races are in
``test_agent_provider_authentication.py`` and
``test_chatgpt_authentication_races.py``. These are the refusals around them:
a credential-vault record that is incomplete or malformed, an OAuth or
metadata response that is unreadable, an identity token that does not verify,
a loopback callback that is not the one this attempt is waiting for, a token
refresh whose answer cannot be used, and an account selection that changes
while it is being checked. Every credential here is synthetic, and no request
leaves the machine: the OpenAI endpoints are replaced at ``_open_https`` and
the browser at ``webbrowser.open``, while the callback is a real loopback
listener driven over HTTP.
"""

from __future__ import annotations

import hashlib
import http.client
import json
import ssl
import sys
import time
import types
import urllib.error
import urllib.parse
from io import BytesIO

import pytest

MODELS_URL = "https://api.openai.com/v1/models"
PLAN_SCOPES = ["openid", "offline_access", "resource.invoke", "chatgpt.tokens.use.direct"]


@pytest.fixture
def chatgpt(monkeypatch, tmp_path):
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
    yield oauth, vault, values
    # A test that leaves an attempt open must not leave its listener running.
    for state in list(oauth._pending):
        try:
            oauth._finish_attempt(state)
        except Exception:  # noqa: BLE001 - teardown only
            oauth._pending.pop(state, None)


def _account(client_id="oaiapp_test", *, expires_at=9_999_999_999, scopes=None,
             subject="subject-test", **extra):
    account = {
        "client_id": client_id, "subject": subject, "issuer": "https://auth.openai.com",
        "email": f"{client_id}@example.invalid", "ext_agent_host_id": "urn:uuid:test-host",
        "id_token": "synthetic-identity", "access_token": f"access-{client_id}",
        "refresh_token": f"refresh-{client_id}", "token_type": "Bearer",
        "expires_at": expires_at, "scopes": list(PLAN_SCOPES if scopes is None else scopes),
        "session_id": f"session-{client_id}",
    }
    account.update(extra)
    return account


def _save(oauth, *accounts, active="oaiapp_test", **extra):
    record = {"host_id": "urn:uuid:test-host", "active_client_id": active,
              "accounts": list(accounts)}
    record.update(extra)
    oauth._write_store(record)


class _Reply(BytesIO):
    def __init__(self, body: bytes, status: int = 200):
        super().__init__(body)
        self.status = status

    def __enter__(self):
        return self

    def __exit__(self, *_args):
        self.close()


def _http_error(url, code, body=b""):
    return urllib.error.HTTPError(url, code, "synthetic", {}, BytesIO(body))


def _raise(error):
    def raising(*_args, **_kwargs):
        raise error
    return raising


def _tokens(**overrides):
    response = {"access_token": "synthetic-access", "refresh_token": "synthetic-refresh",
                "id_token": "synthetic-identity", "token_type": "Bearer",
                "expires_in": 3600, "scope": " ".join(PLAN_SCOPES)}
    response.update(overrides)
    return response


def _begin(oauth, monkeypatch, **kwargs):
    """Start an attempt with a stand-in browser; return its state and listener port."""
    opened = []
    monkeypatch.setattr(oauth.webbrowser, "open",
                        lambda url, **_kw: opened.append(url) or True)
    result = oauth.connect(**kwargs)
    assert result["ok"] is True, result
    query = urllib.parse.parse_qs(urllib.parse.urlsplit(opened[-1]).query)
    port = urllib.parse.urlsplit(query["redirect_uri"][0]).port
    return query["state"][0], port, query


def _get(port, path):
    connection = http.client.HTTPConnection("127.0.0.1", port, timeout=10)
    try:
        connection.request("GET", path)
        response = connection.getresponse()
        return response.status, response.read().decode("utf-8")
    finally:
        connection.close()


def _callback(port, **params):
    return _get(port, "/auth/callback?" + urllib.parse.urlencode(params))


# -- the credential vault --------------------------------------------------


def test_no_secure_vault_means_no_sign_in_and_nothing_is_written(monkeypatch, chatgpt):
    from fastmdxplora.agent import models

    oauth, _vault, values = chatgpt
    monkeypatch.setattr(models, "_secure_keyring", lambda: None)
    assert oauth._read_store() == {"host_id": "", "active_client_id": "", "accounts": []}
    with pytest.raises(oauth.CredentialStoreError, match="No supported OS credential vault"):
        oauth._write_store({"accounts": []})
    answer = oauth.status()
    assert answer["available"] is False and answer["connected"] is False
    assert answer["accounts"] == [] and "API key" in answer["message"]
    started = oauth.connect()
    assert started["ok"] is False and "credential vault is required" in started["message"]
    assert values == {} and oauth._pending == {}


@pytest.mark.parametrize("manifest", [
    {"generation": "not-a-hex-generation", "count": 1},
    {"generation": "0" * 31 + "g", "count": 1},
    {"generation": "0" * 32, "count": 0},
    {"generation": "0" * 32, "count": 4097},
    {"generation": "0" * 32, "count": "2"},
])
def test_a_malformed_vault_manifest_is_refused_without_exposing_it(chatgpt, manifest):
    oauth, vault, _values = chatgpt
    raw = json.dumps({"vault_chunks": manifest})
    with pytest.raises(oauth.CredentialStoreError, match="manifest is invalid"):
        oauth._chunk_keys(raw)
    vault.set_password(oauth._SERVICE, oauth._vault_user(), raw)
    with pytest.raises(oauth.CredentialStoreError, match="could not read") as caught:
        oauth._read_store()
    assert "not-a-hex" not in str(caught.value)
    assert caught.value.__cause__ is None and caught.value.__suppress_context__


def test_a_vault_record_missing_a_part_is_unreadable_not_empty(chatgpt):
    oauth, vault, values = chatgpt
    _save(oauth, _account())
    part = next(key for key in values if key[0] != oauth._SERVICE)
    vault.delete_password(*part)
    with pytest.raises(oauth.CredentialStoreError, match="could not read"):
        oauth._read_store()


def test_a_vault_record_that_is_not_an_object_reads_as_no_sign_in(chatgpt):
    oauth, vault, _values = chatgpt
    vault.set_password(oauth._SERVICE, oauth._vault_user(), json.dumps([1, 2, 3]))
    store = oauth._read_store()
    assert store["accounts"] == [] and store["active_client_id"] == ""
    assert store["auth_epoch"] == 0 and store["pending_operation"] == ""


def test_a_failed_rollback_still_reports_the_failed_write_and_keeps_the_old_record(chatgpt):
    oauth, vault, values = chatgpt
    _save(oauth, _account())
    before = dict(values)
    set_password, delete_password = vault.set_password, vault.delete_password
    deleted = []

    def refuse_manifest(service, username, value):
        if service == oauth._SERVICE:
            raise OSError("synthetic vault refusal")
        set_password(service, username, value)

    def refuse_delete(service, username):
        deleted.append(service)
        raise OSError("synthetic delete refusal")

    vault.set_password, vault.delete_password = refuse_manifest, refuse_delete
    with pytest.raises(oauth.CredentialStoreError, match="could not save") as caught:
        oauth._write_store({"active_client_id": "replacement", "accounts": []})
    assert "synthetic" not in str(caught.value)
    # Every new part was offered for removal even though the vault refused it.
    assert deleted and all(name.startswith(oauth._SERVICE + ":") for name in deleted)
    vault.set_password, vault.delete_password = set_password, delete_password
    assert values[(oauth._SERVICE, oauth._vault_user())] == before[
        (oauth._SERVICE, oauth._vault_user())]
    assert oauth._read_store()["active_client_id"] == "oaiapp_test"


def test_windows_verifies_openai_with_the_system_trust_store(monkeypatch):
    from fastmdxplora.agent import openai_chatgpt as oauth

    made = []
    truststore = types.ModuleType("truststore")
    truststore.SSLContext = lambda protocol: made.append(protocol) or "system-trust-context"
    monkeypatch.setitem(sys.modules, "truststore", truststore)
    monkeypatch.setattr(oauth, "sys", types.SimpleNamespace(platform="win32"))
    assert oauth._https_context() == "system-trust-context"
    assert made == [ssl.PROTOCOL_TLS_CLIENT]


def test_windows_process_lock_locks_and_unlocks_one_byte(monkeypatch, chatgpt):
    oauth, _vault, _values = chatgpt
    calls = []
    msvcrt = types.ModuleType("msvcrt")
    msvcrt.LK_LOCK, msvcrt.LK_UNLCK = 1, 0
    msvcrt.locking = lambda fileno, mode, size: calls.append((mode, size))
    monkeypatch.setitem(sys.modules, "msvcrt", msvcrt)
    monkeypatch.setattr(oauth, "os", types.SimpleNamespace(name="nt"))
    with oauth._process_lock(purpose="store"):
        assert calls == [(msvcrt.LK_LOCK, 1)]
    assert calls == [(msvcrt.LK_LOCK, 1), (msvcrt.LK_UNLCK, 1)]


# -- HTTPS responses from OpenAI ---------------------------------------------


@pytest.mark.parametrize("url", [
    "http://auth.openai.com/api/accounts/oauth/token",
    "https://auth.openai.com.example.invalid/token",
    "https://user:secret@auth.openai.com/token",
    "https://auth.openai.com:8443/token",
])
def test_credentials_are_sent_only_to_the_documented_https_hosts(url):
    import urllib.request

    from fastmdxplora.agent import openai_chatgpt as oauth

    with pytest.raises(oauth.SignInError, match="unexpected HTTPS endpoint") as caught:
        oauth._open_https(urllib.request.Request(url), timeout=1)
    assert "secret" not in str(caught.value)


def test_transport_failures_are_described_without_their_details():
    from fastmdxplora.agent import openai_chatgpt as oauth

    assert "timed out" in oauth._transport_message(TimeoutError("synthetic"))
    assert "timed out" in oauth._transport_message(urllib.error.URLError(TimeoutError()))
    refused = oauth._transport_message(urllib.error.URLError(ConnectionRefusedError("x")))
    assert "could not be reached" in refused and "x" not in refused.split()


def test_a_token_response_is_returned_only_when_it_is_a_json_object(monkeypatch):
    from fastmdxplora.agent import openai_chatgpt as oauth

    sent = []

    def reply(body):
        def opened(request, **kwargs):
            sent.append((request, kwargs))
            return _Reply(body)
        return opened

    monkeypatch.setattr(oauth, "_open_https", reply(json.dumps({"access_token": "a"}).encode()))
    assert oauth._exchange({"grant_type": "refresh_token", "refresh_token": "r"},
                           timeout=7) == {"access_token": "a"}
    request, kwargs = sent[0]
    assert request.full_url == oauth.TOKEN_URL and request.get_method() == "POST"
    assert urllib.parse.parse_qs(request.data.decode()) == {
        "grant_type": ["refresh_token"], "refresh_token": ["r"]}
    assert kwargs == {"timeout": 7}
    for body in (b"\xff\xfe", b"not json", b"[1, 2]", b'"text"'):
        monkeypatch.setattr(oauth, "_open_https", reply(body))
        with pytest.raises(oauth.SignInError, match="unreadable response"):
            oauth._exchange({"grant_type": "refresh_token"})


@pytest.mark.parametrize("body,code", [
    (b"<html>synthetic-secret</html>", 500),
    (json.dumps({"error": "invalid_request", "echo": "synthetic-secret"}).encode(), 400),
    (json.dumps({"error": {"code": 7}}).encode(), 400),
    (json.dumps(["synthetic-secret"]).encode(), 400),
])
def test_a_rejected_token_request_reports_only_its_status(monkeypatch, body, code):
    from fastmdxplora.agent import openai_chatgpt as oauth

    monkeypatch.setattr(oauth, "_open_https",
                        _raise(_http_error(oauth.TOKEN_URL, code, body)))
    with pytest.raises(oauth.SignInError) as caught:
        oauth._exchange({"code": "synthetic-secret"})
    assert not isinstance(caught.value, oauth.InvalidGrant)
    assert str(caught.value) == f"OpenAI sign-in returned HTTP {code}. Try again."


def test_an_unreachable_token_endpoint_is_a_sign_in_error(monkeypatch):
    from fastmdxplora.agent import openai_chatgpt as oauth

    monkeypatch.setattr(oauth, "_open_https", _raise(ConnectionResetError("synthetic")))
    with pytest.raises(oauth.SignInError, match="could not be reached"):
        oauth._exchange({"grant_type": "authorization_code"})


def test_identity_metadata_must_be_a_json_object(monkeypatch):
    from fastmdxplora.agent import openai_chatgpt as oauth

    monkeypatch.setattr(oauth, "_open_https",
                        lambda request, **_kw: _Reply(b'{"issuer": "https://auth.openai.com"}'))
    assert oauth._fetch_json(oauth.DISCOVERY_URL) == {"issuer": oauth.ISSUER}
    for body in (b"\xff", b"{broken", b"[]"):
        monkeypatch.setattr(oauth, "_open_https", lambda request, _b=body, **_kw: _Reply(_b))
        with pytest.raises(oauth.SignInError, match="metadata was not valid JSON"):
            oauth._fetch_json(oauth.DISCOVERY_URL)
    monkeypatch.setattr(oauth, "_open_https", _raise(TimeoutError()))
    with pytest.raises(oauth.SignInError, match="timed out"):
        oauth._fetch_json(oauth.DISCOVERY_URL)


# -- the identity token ------------------------------------------------------

_MODULUS = int(
    "9d0b9558969c6516a0614e724784cd45aab568f23cfc7f408f0f02fbbefee9da"
    "9d6ba0aa067f6d1ab0e28699f77810791b4840edb80614873b532c359f98b480"
    "6b9aac9b12a94150300b393fd4a872f6703d19e0fd09faa7fa1e528b52b241d"
    "ab4171ed17219a208fa0178652917c583708d344e248eb1c49233016c3ab66503", 16)
_PRIVATE = int(
    "419616d71a7f525072c670a6956bc6f6b36673dfc6cc84b0878e736683cb9e0a"
    "a8f7cbef762b3d7b8a27f0a2d891120021feb8ef0e8dab8116ab93a9097012e"
    "16db7393e6f0acbb6479fded6ac595ae6a9704da188f9d8cafab4a5cd1c1620"
    "c65a968435aeeb1fceac4f249ef26aa7d59afeb999fff2b3102b0bbb9cb4f40941", 16)
_KID = "synthetic-signing-key"


def _sign(oauth, message: str) -> str:
    digest = hashlib.sha256(message.encode("ascii")).digest()
    info = bytes.fromhex("3031300d060960864801650304020105000420") + digest
    size = (_MODULUS.bit_length() + 7) // 8
    block = b"\x00\x01" + b"\xff" * (size - len(info) - 3) + b"\x00" + info
    return oauth._b64url(pow(int.from_bytes(block, "big"), _PRIVATE, _MODULUS)
                         .to_bytes(size, "big"))


def _identity(oauth, *, header=None, claims=None):
    header = header or {"alg": "RS256", "kid": _KID}
    claims = claims or {"iss": oauth.ISSUER, "sub": "synthetic-subject", "aud": "oaiapp_test",
                        "exp": int(time.time()) + 3600, "nonce": "synthetic-nonce"}
    encoded = ".".join(oauth._b64url(json.dumps(part).encode()) for part in (header, claims))
    return encoded + "." + _sign(oauth, encoded)


def _jwk(oauth, **overrides):
    key = {"kid": _KID, "kty": "RSA", "use": "sig", "alg": "RS256",
           "n": oauth._b64url(_MODULUS.to_bytes((_MODULUS.bit_length() + 7) // 8, "big")),
           "e": oauth._b64url((65537).to_bytes(3, "big"))}
    key.update(overrides)
    return {key_name: value for key_name, value in key.items() if value is not None}


def _metadata(monkeypatch, oauth, *, discovery=None, keys=None):
    discovery = discovery or {"issuer": oauth.ISSUER,
                              "jwks_uri": oauth.ISSUER + "/.well-known/jwks.json"}
    keys = {"keys": [_jwk(oauth)]} if keys is None else keys
    asked = []

    def fetch(url, **_kwargs):
        asked.append(url)
        return discovery if url == oauth.DISCOVERY_URL else keys

    monkeypatch.setattr(oauth, "_fetch_json", fetch)
    return asked


def _verify(oauth, token):
    return oauth._verify_id_token(token, client_id="oaiapp_test", nonce="synthetic-nonce")


def test_identity_tokens_that_are_not_three_rs256_parts_are_refused_before_any_request(
        monkeypatch):
    from fastmdxplora.agent import openai_chatgpt as oauth

    asked = _metadata(monkeypatch, oauth)
    with pytest.raises(oauth.SignInError, match="invalid identity token"):
        _verify(oauth, "only.two")
    with pytest.raises(oauth.SignInError, match="unsupported identity-token signature"):
        _verify(oauth, _identity(oauth, header={"alg": "HS256", "kid": _KID}))
    assert asked == []
    assert _verify(oauth, _identity(oauth))["sub"] == "synthetic-subject"


@pytest.mark.parametrize("discovery,keys,message", [
    ({"issuer": "https://issuer.example.invalid",
      "jwks_uri": "https://auth.openai.com/jwks"}, None, "unexpected issuer"),
    ({"issuer": "https://auth.openai.com",
      "jwks_uri": "https://keys.example.invalid/jwks"}, None, "unexpected signing-key URL"),
    (None, {"keys": "not-a-list"}, "signing key was not recognized"),
    (None, {"keys": [{"kid": "another-key", "kty": "RSA"}]}, "signing key was not recognized"),
    (None, {"keys": [{"kid": _KID, "kty": "EC"}]}, "signing key was not recognized"),
    (None, {"keys": [{"kid": _KID, "kty": "RSA", "e": "AQAB"}]}, "signing key was malformed"),
])
def test_identity_metadata_and_signing_keys_must_be_openai_s(monkeypatch, discovery, keys,
                                                              message):
    from fastmdxplora.agent import openai_chatgpt as oauth

    _metadata(monkeypatch, oauth, discovery=discovery, keys=keys)
    with pytest.raises(oauth.SignInError, match=message):
        _verify(oauth, _identity(oauth))


def test_a_full_size_signature_over_different_claims_does_not_validate(monkeypatch):
    from fastmdxplora.agent import openai_chatgpt as oauth

    _metadata(monkeypatch, oauth)
    token = _identity(oauth)
    header, claims, _signature = token.split(".")
    forged = oauth._b64url(json.dumps({"iss": oauth.ISSUER, "sub": "someone-else",
                                       "aud": "oaiapp_test", "nonce": "synthetic-nonce",
                                       "exp": int(time.time()) + 3600}).encode())
    with pytest.raises(oauth.SignInError, match="signature did not validate"):
        _verify(oauth, f"{header}.{forged}.{_sign(oauth, header + '.' + claims)}")


def test_identity_claims_without_a_subject_or_with_a_future_start_are_refused(monkeypatch):
    from fastmdxplora.agent import openai_chatgpt as oauth

    _metadata(monkeypatch, oauth)
    now = int(time.time())
    base = {"iss": oauth.ISSUER, "aud": ["other", "oaiapp_test"], "exp": now + 600,
            "nonce": "synthetic-nonce"}
    assert _verify(oauth, _identity(oauth, claims={**base, "sub": "listed-audience"}))[
        "sub"] == "listed-audience"
    for claims in ({**base}, {**base, "sub": "s", "nbf": now + 3600},
                   {**base, "sub": "s", "aud": 7}):
        with pytest.raises(oauth.SignInError, match="did not validate"):
            _verify(oauth, _identity(oauth, claims=claims))


@pytest.mark.parametrize("response,message", [
    (_tokens(refresh_token=""), "tokens needed for a renewable sign-in"),
    (_tokens(id_token=None), "tokens needed for a renewable sign-in"),
    (_tokens(expires_in="later"), "valid token expiry"),
    (_tokens(expires_in=0), "valid token expiry"),
    (_tokens(token_type="MAC"), "unsupported token type"),
])
def test_token_responses_must_be_renewable_bearer_tokens_with_an_expiry(
        monkeypatch, response, message):
    from fastmdxplora.agent import openai_chatgpt as oauth

    monkeypatch.setattr(oauth, "_verify_id_token",
                        lambda *_a, **_kw: {"sub": "synthetic-subject"})
    with pytest.raises(oauth.SignInError, match=message):
        oauth._validate_tokens(response, client_id="oaiapp_test", nonce="n")
    accepted = oauth._validate_tokens(_tokens(scope="b a a"), client_id="oaiapp_test",
                                      nonce="n")
    assert accepted["scopes"] == ["a", "b"] and accepted["email"] == ""


# -- the loopback callback ---------------------------------------------------


def test_requests_that_are_not_this_attempt_s_callback_leave_it_waiting(
        monkeypatch, chatgpt):
    oauth, _vault, _values = chatgpt
    exchanged = []
    monkeypatch.setattr(oauth, "_exchange", lambda form: exchanged.append(form) or {})
    state, port, _query = _begin(oauth, monkeypatch)
    status, body = _get(port, "/favicon.ico")
    assert status == 404 and "not the sign-in callback" in body
    status, body = _get(port, "/auth/callback?state=a&state=b&code=c")
    assert status == 400 and "repeated parameters" in body
    status, body = _callback(port, state="an-unknown-state", code="c")
    assert status == 400 and "expired or did not match" in body
    # An attempt registered by a different listener is not this one's to finish.
    foreign = {"server": object(), "processing": False}
    oauth._pending["foreign-state"] = foreign
    status, body = _callback(port, state="foreign-state", code="c")
    assert status == 400 and "expired or did not match" in body
    assert foreign["processing"] is False
    del oauth._pending["foreign-state"]
    assert exchanged == [] and state in oauth._pending
    assert "history.replaceState" not in _get(port, "/elsewhere")[1]


@pytest.mark.parametrize("params,message", [
    ({"error": "server_error"}, "was not completed"),
    ({"code": ""}, "did not include an authorization code"),
    ({"code": "c"}, "did not return the newly registered client ID"),
    ({"code": "c", "client_id": "dynamic_agent_client"},
     "did not return the newly registered client ID"),
])
def test_a_callback_without_a_usable_code_or_client_ends_the_attempt(
        monkeypatch, chatgpt, params, message):
    oauth, _vault, _values = chatgpt
    exchanged = []
    monkeypatch.setattr(oauth, "_exchange", lambda form: exchanged.append(form) or {})
    state, port, _query = _begin(oauth, monkeypatch)
    status, body = _callback(port, state=state, **params)
    assert status == 400 and message in body
    assert exchanged == [] and oauth._pending == {}
    assert oauth.status()["message"] == oauth._pending_note
    assert message in oauth._pending_note
    store = oauth._read_store()
    assert store["pending_operation"] == "" and store["accounts"] == []


def test_a_callback_for_another_registration_is_refused_for_a_saved_account(
        monkeypatch, chatgpt):
    oauth, _vault, _values = chatgpt
    _save(oauth, _account(scopes=["openid"]))
    exchanged = []
    monkeypatch.setattr(oauth, "_exchange", lambda form: exchanged.append(form) or {})
    state, port, query = _begin(oauth, monkeypatch, client_id="oaiapp_test")
    assert query["client_id"] == ["oaiapp_test"]
    assert query["login_hint"] == ["oaiapp_test@example.invalid"]
    assert query["id_token_hint"] == ["synthetic-identity"]
    assert query["prompt"] == ["consent"]
    status, body = _callback(port, state=state, code="c", client_id="oaiapp_other")
    assert status == 400 and "different account registration" in body
    assert exchanged == []
    assert oauth._read_store()["accounts"][0]["access_token"] == "access-oaiapp_test"


def test_a_saved_client_id_cannot_be_taken_over_by_another_chatgpt_account(
        monkeypatch, chatgpt):
    oauth, _vault, _values = chatgpt
    _save(oauth, _account())
    monkeypatch.setattr(oauth, "_exchange", lambda _form: _tokens())
    monkeypatch.setattr(oauth, "_verify_id_token", lambda *_a, **_kw: {
        "sub": "a-different-subject", "email": "other@example.invalid"})
    state, port, query = _begin(oauth, monkeypatch, client_id="oaiapp_test")
    assert "prompt" not in query
    status, body = _callback(port, state=state, code="c")
    assert status == 400 and "belongs to another ChatGPT account" in body
    saved = oauth._read_store()["accounts"]
    assert len(saved) == 1 and saved[0]["subject"] == "subject-test"
    assert saved[0]["access_token"] == "access-oaiapp_test"


def test_a_callback_whose_attempt_was_replaced_during_the_exchange_saves_nothing(
        monkeypatch, chatgpt):
    oauth, _vault, _values = chatgpt
    monkeypatch.setattr(oauth, "_verify_id_token", lambda *_a, **_kw: {"sub": "s"})

    def exchange(_form):
        # Another tab replaces this intent while the token request is out.
        oauth._pending[state] = dict(oauth._pending[state])
        return _tokens()

    monkeypatch.setattr(oauth, "_exchange", exchange)
    state, port, _query = _begin(oauth, monkeypatch)
    status, body = _callback(port, state=state, code="c", client_id="oaiapp_new")
    assert status == 400 and "was canceled" in body
    store = oauth._read_store()
    assert store["accounts"] == [] and store["active_client_id"] == ""
    assert store["pending_client_id"] == "oaiapp_new"


def test_an_unexpected_failure_while_validating_is_reported_without_its_detail(
        monkeypatch, chatgpt):
    oauth, _vault, _values = chatgpt
    monkeypatch.setattr(oauth, "_exchange", _raise(KeyError("synthetic-secret-detail")))
    state, port, _query = _begin(oauth, monkeypatch)
    status, body = _callback(port, state=state, code="c", client_id="oaiapp_new")
    assert status == 400 and "failed while validating the response" in body
    assert "synthetic-secret-detail" not in body + oauth._pending_note
    assert oauth._pending == {} and oauth._read_store()["accounts"] == []


def test_a_vault_failure_while_ending_the_attempt_is_what_the_browser_is_told(
        monkeypatch, chatgpt):
    oauth, _vault, _values = chatgpt
    monkeypatch.setattr(oauth, "_exchange", lambda _form: _tokens())
    monkeypatch.setattr(oauth, "_verify_id_token", lambda *_a, **_kw: {"sub": "s"})
    finish = oauth._finish_attempt

    def finish_then_fail(state):
        finish(state)
        raise oauth.CredentialStoreError("The vault locked while finishing sign-in.")

    monkeypatch.setattr(oauth, "_finish_attempt", finish_then_fail)
    state, port, _query = _begin(oauth, monkeypatch)
    status, body = _callback(port, state=state, code="c", client_id="oaiapp_new")
    assert status == 400 and "vault locked while finishing" in body
    assert oauth._pending_note == "The vault locked while finishing sign-in."
    assert oauth._pending == {}


# -- starting and abandoning an attempt -------------------------------------


def test_only_one_attempt_waits_and_only_a_saved_account_can_be_selected(
        monkeypatch, chatgpt):
    oauth, _vault, _values = chatgpt
    unknown = oauth.connect(client_id="oaiapp_unknown")
    assert unknown == {"ok": False,
                       "message": "Select a saved ChatGPT account, then try again."}
    state, _port, _query = _begin(oauth, monkeypatch)
    second = oauth.connect()
    assert second["ok"] is False and "already waiting" in second["message"]
    assert list(oauth._pending) == [state]
    answer = oauth.status()
    assert answer["pending"] is True and answer["connected"] is None
    assert "Waiting for ChatGPT sign-in" in answer["message"]


def test_a_callback_listener_that_cannot_start_leaves_no_attempt(monkeypatch, chatgpt):
    oauth, _vault, values = chatgpt
    monkeypatch.setattr(oauth, "_LoopbackServer", _raise(OSError("address in use")))
    answer = oauth.connect()
    assert answer["ok"] is False and "callback could not be started" in answer["message"]
    assert oauth._pending == {} and values == {}


def test_an_attempt_superseded_while_starting_is_abandoned(monkeypatch, chatgpt):
    oauth, _vault, _values = chatgpt
    closed = []

    class Superseded(oauth._LoopbackServer):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            store = oauth._read_store()
            store["auth_epoch"] += 1  # another dashboard started first
            oauth._write_store(store)

        def server_close(self):
            closed.append(True)
            super().server_close()

    monkeypatch.setattr(oauth, "_LoopbackServer", Superseded)
    answer = oauth.connect()
    assert answer == {"ok": False,
                      "message": "ChatGPT sign-in changed while starting. Try again."}
    assert closed == [True] and oauth._pending == {}
    assert oauth._read_store()["pending_operation"] == ""


def test_an_attempt_canceled_before_it_is_registered_is_not_opened(monkeypatch, chatgpt):
    oauth, _vault, _values = chatgpt
    write = oauth._write_store
    opened = []

    def write_then_cancel(record):
        write(record)
        if record.get("pending_operation"):
            # A disconnect in another process lands between the two steps.
            monkeypatch.setattr(oauth, "_write_store", write)
            write({**record, "pending_operation": "", "auth_epoch": record["auth_epoch"] + 1})

    monkeypatch.setattr(oauth, "_write_store", write_then_cancel)
    monkeypatch.setattr(oauth.webbrowser, "open", lambda url, **_kw: opened.append(url))
    answer = oauth.connect()
    assert answer["ok"] is False and "canceled or superseded" in answer["message"]
    assert opened == [] and oauth._pending == {}


@pytest.mark.parametrize("browser", [lambda *_a, **_kw: False,
                                     _raise(RuntimeError("no display"))])
def test_a_browser_that_does_not_open_ends_the_attempt(monkeypatch, chatgpt, browser):
    oauth, _vault, _values = chatgpt
    monkeypatch.setattr(oauth.webbrowser, "open", browser)
    answer = oauth.connect()
    assert answer["ok"] is False and "system browser could not open" in answer["message"]
    assert oauth._pending == {}
    store = oauth._read_store()
    assert store["pending_operation"] == "" and store["auth_epoch"] == 2


def test_an_abandoned_attempt_expires_once_and_says_so(monkeypatch, chatgpt):
    oauth, _vault, _values = chatgpt
    state, port, _query = _begin(oauth, monkeypatch)
    oauth._expire_attempt(state)
    note = "ChatGPT sign-in timed out or was closed. Start it again when ready."
    assert oauth._pending == {} and oauth._pending_note == note
    assert oauth.status()["message"] == note
    oauth._pending_note = "a later note"
    oauth._expire_attempt(state)
    assert oauth._pending_note == "a later note"
    deadline = time.monotonic() + 5
    while time.monotonic() < deadline:
        try:
            _get(port, "/auth/callback")
        except OSError:
            break
        time.sleep(0.05)
    else:
        pytest.fail("the expired attempt's listener kept accepting connections")


# -- refreshing, status, models, disconnecting and switching -----------------


@pytest.mark.parametrize("account,message", [
    (_account(scopes=["openid", "offline_access"]), "plan use is not authorized"),
    (_account(token_type="MAC"), "does not use bearer tokens"),
    (_account(refresh_token=""), "No active ChatGPT sign-in"),
])
def test_a_saved_session_is_used_only_as_a_plan_bearer_session(chatgpt, account, message):
    oauth, _vault, _values = chatgpt
    _save(oauth, account)
    with pytest.raises(oauth.SignInError, match=message):
        oauth.access_token()


def test_a_token_refreshed_by_another_process_is_used_without_refreshing_again(
        monkeypatch, chatgpt):
    oauth, _vault, _values = chatgpt
    _save(oauth, _account(expires_at=0))
    read = oauth._read_store
    reads = []

    def read_after_other_refresh():
        store = read()
        reads.append(True)
        if len(reads) == 2:
            store["accounts"][0].update(access_token="refreshed-elsewhere",
                                        expires_at=int(time.time()) + 3600)
        return store

    monkeypatch.setattr(oauth, "_read_store", read_after_other_refresh)
    monkeypatch.setattr(oauth, "_exchange", _raise(AssertionError("must not refresh")))
    assert oauth.access_token() == "refreshed-elsewhere"


def test_a_refresh_without_an_access_token_keeps_the_saved_session(monkeypatch, chatgpt):
    oauth, _vault, _values = chatgpt
    _save(oauth, _account(expires_at=0))
    monkeypatch.setattr(oauth, "_exchange", lambda _form: {"refresh_token": "rotated"})
    with pytest.raises(oauth.SignInError, match="returned no access token"):
        oauth.access_token()
    saved = oauth._read_store()
    assert saved["active_client_id"] == "oaiapp_test"
    assert saved["accounts"][0]["refresh_token"] == "refresh-oaiapp_test"


@pytest.mark.parametrize("response", [
    {"access_token": "new", "expires_in": "soon"},
    {"access_token": "new", "expires_in": -5},
    {"access_token": "new", "expires_in": 3600, "token_type": "MAC"},
])
def test_unusable_refreshed_credentials_are_cleared_and_deselected(
        monkeypatch, chatgpt, response):
    oauth, _vault, _values = chatgpt
    _save(oauth, _account(expires_at=0))
    monkeypatch.setattr(oauth, "_exchange", lambda _form: dict(response))
    with pytest.raises(oauth.SignInError, match="invalid refreshed credentials"):
        oauth.access_token()
    saved = oauth._read_store()
    assert saved["active_client_id"] == ""
    assert saved["accounts"][0]["client_id"] == "oaiapp_test"
    assert not any(key in saved["accounts"][0]
                   for key in ("access_token", "refresh_token", "scopes"))


def test_a_refresh_that_drops_plan_permission_keeps_the_rotated_token_but_deselects(
        monkeypatch, chatgpt):
    oauth, _vault, _values = chatgpt
    _save(oauth, _account(expires_at=0))
    monkeypatch.setattr(oauth, "_exchange", lambda _form: {
        "access_token": "rotated-access", "refresh_token": "rotated-refresh",
        "expires_in": 3600, "scope": "openid offline_access"})
    with pytest.raises(oauth.SignInError, match="no longer includes ChatGPT plan permission"):
        oauth.access_token()
    saved = oauth._read_store()
    assert saved["active_client_id"] == ""
    account = saved["accounts"][0]
    assert account["refresh_token"] == "rotated-refresh"
    assert account["scopes"] == ["offline_access", "openid"]
    assert oauth.account_choices(saved)[0]["connected"] is False


def test_status_reports_a_failed_check_and_an_account_without_plan_use(monkeypatch, chatgpt):
    oauth, _vault, _values = chatgpt
    _save(oauth, _account())
    monkeypatch.setattr(oauth, "access_token",
                        _raise(oauth.SignInError("OpenAI sign-in could not be reached.")))
    answer = oauth.status()
    assert answer["connected"] is False
    assert answer["message"] == "OpenAI sign-in could not be reached."
    assert answer["accounts"][0]["active"] is True
    _save(oauth, _account(scopes=["openid"]))
    answer = oauth.status()
    assert answer["connected"] is False
    assert "plan use was not authorized" in answer["message"]
    assert "oaiapp_test@example.invalid" in answer["message"]


def test_status_does_not_confirm_an_account_switched_while_it_was_checked(
        monkeypatch, chatgpt):
    oauth, _vault, _values = chatgpt
    _save(oauth, _account(), _account("oaiapp_other"))

    def switch_during_check():
        store = oauth._read_store()
        store["active_client_id"] = "oaiapp_other"
        oauth._write_store(store)
        return "access-oaiapp_test"

    monkeypatch.setattr(oauth, "access_token", switch_during_check)
    answer = oauth.status()
    assert answer["connected"] is None
    assert answer["message"] == ("The ChatGPT account changed while checking it. "
                                 "Refresh its status.")
    assert answer["active_client_id"] == "oaiapp_other"


@pytest.mark.parametrize("error,message", [
    (_http_error(MODELS_URL, 401), "rejected the active"),
    (_http_error(MODELS_URL, 403), "policy or account permissions"),
    (_http_error(MODELS_URL, 503), "temporarily unavailable"),
    (_http_error(MODELS_URL, 418), "returned HTTP 418"),
    (urllib.error.URLError("synthetic"), "could not be reached"),
    (TimeoutError(), "could not be reached"),
])
def test_model_listing_failures_say_what_to_do(monkeypatch, error, message):
    from fastmdxplora.agent import openai_chatgpt as oauth

    monkeypatch.setattr(oauth, "access_token", lambda: "synthetic-access")
    monkeypatch.setattr(oauth, "_open_https", _raise(error))
    with pytest.raises(oauth.SignInError, match=message) as caught:
        oauth.list_models()
    assert "synthetic-access" not in str(caught.value)


@pytest.mark.parametrize("body", [b'{"data": []}', b'["gpt"]', b'{"models": "gpt"}'])
def test_an_unrecognized_model_catalog_is_refused(monkeypatch, body):
    from fastmdxplora.agent import openai_chatgpt as oauth

    monkeypatch.setattr(oauth, "access_token", lambda: "synthetic-access")
    monkeypatch.setattr(oauth, "_open_https", lambda _request, **_kw: _Reply(body))
    with pytest.raises(oauth.SignInError, match="unrecognized ChatGPT model catalog"):
        oauth.list_models()


def test_an_unreachable_revocation_endpoint_is_retried_then_reported(monkeypatch, chatgpt):
    oauth, _vault, _values = chatgpt
    _save(oauth, _account())
    monkeypatch.setattr(oauth, "_fetch_json", lambda _url: {
        "revocation_endpoint": oauth.ISSUER + "/api/accounts/oauth/revoke"})
    waits, attempts = [], []
    monkeypatch.setattr(oauth.time, "sleep", waits.append)

    def unreachable(request, **_kwargs):
        attempts.append(request.full_url)
        raise urllib.error.URLError(ConnectionResetError("synthetic"))

    monkeypatch.setattr(oauth, "_open_https", unreachable)
    answer = oauth.disconnect()
    assert answer["ok"] is True and answer["remote_revoked"] is None
    assert "could not be reached for remote sign-out" in answer["message"]
    assert len(attempts) == 3 and waits == [1, 2]
    assert "refresh_token" not in oauth._read_store()["accounts"][0]


def test_revocation_metadata_outside_openai_is_not_used(monkeypatch, chatgpt):
    oauth, _vault, _values = chatgpt
    _save(oauth, _account())
    monkeypatch.setattr(oauth, "_fetch_json", lambda _url: {
        "revocation_endpoint": "https://revoke.example.invalid/"})
    monkeypatch.setattr(oauth, "_open_https", _raise(AssertionError("must not be sent")))
    answer = oauth.disconnect()
    assert answer["remote_revoked"] is None
    assert answer["message"].endswith("OpenAI revocation metadata was not valid.")


def test_disconnecting_an_account_with_no_renewable_session_skips_revocation(
        monkeypatch, chatgpt):
    oauth, _vault, _values = chatgpt
    _save(oauth, _account(refresh_token=""), active="")
    monkeypatch.setattr(oauth, "_fetch_json", _raise(AssertionError("no revocation")))
    answer = oauth.disconnect(client_id="oaiapp_test")
    assert answer["ok"] is True and answer["remote_revoked"] is None
    assert answer["message"].endswith(
        "No local renewable session was available to revoke remotely.")


def test_only_a_saved_account_with_a_plan_grant_can_be_activated(chatgpt):
    oauth, _vault, _values = chatgpt
    _save(oauth, _account(), _account("oaiapp_noplan", scopes=["openid"]),
          _account("oaiapp_signedout", refresh_token=""))
    assert oauth.activate(client_id="oaiapp_missing") == {
        "ok": False, "message": "Select a saved ChatGPT account."}
    for client_id in ("oaiapp_noplan", "oaiapp_signedout"):
        answer = oauth.activate(client_id=client_id)
        assert answer["ok"] is False and "no active ChatGPT plan grant" in answer["message"]
    assert oauth._read_store()["active_client_id"] == "oaiapp_test"


def test_activating_an_account_ends_a_waiting_sign_in(monkeypatch, chatgpt):
    oauth, _vault, _values = chatgpt
    _save(oauth, _account(), _account("oaiapp_other"))
    _state, _port, _query = _begin(oauth, monkeypatch)
    assert oauth.activate(client_id="oaiapp_other") == {
        "ok": True, "message": "The selected ChatGPT account is active."}
    assert oauth._pending == {}
    store = oauth._read_store()
    assert store["active_client_id"] == "oaiapp_other" and store["pending_operation"] == ""


def test_a_newer_selection_during_activation_is_not_overridden(monkeypatch, chatgpt):
    oauth, _vault, _values = chatgpt
    _save(oauth, _account(), _account("oaiapp_other"))

    def newer_selection():
        store = oauth._read_store()
        store.update(active_client_id="oaiapp_test", auth_epoch=store["auth_epoch"] + 1)
        oauth._write_store(store)
        return "access-oaiapp_other"

    monkeypatch.setattr(oauth, "access_token", newer_selection)
    assert oauth.activate(client_id="oaiapp_other") == {
        "ok": False, "message": "ChatGPT account selection changed. Retry."}
    assert oauth._read_store()["active_client_id"] == "oaiapp_test"


def test_a_responses_message_item_is_text_not_a_tool_call_and_cached_input_is_counted(
        monkeypatch):
    from fastmdxplora.agent import turns
    from fastmdxplora.agent import openai_chatgpt as oauth

    events = [
        {"type": "response.output_item.added",
         "item": {"id": "msg_1", "type": "message", "role": "assistant", "content": []}},
        {"type": "response.output_text.delta", "delta": "Frame 7 has the lowest RMSD."},
        {"type": "response.output_item.done",
         "item": {"id": "msg_1", "type": "message", "role": "assistant", "content": []}},
        {"type": "response.completed", "response": {
            "status": "completed",
            "output": [{"id": "msg_1", "type": "message"}, "not-an-item"],
            "usage": {"input_tokens": 30, "output_tokens": 6,
                      "input_tokens_details": {"cached_tokens": 20}}}},
    ]

    class Reply:
        def __enter__(self):
            return self

        def __exit__(self, *_args):
            return None

        def __iter__(self):
            return iter([("data: " + json.dumps(event) + "\n\n").encode()
                         for event in events])

    monkeypatch.setattr(oauth, "_open_https", lambda *_a, **_kw: Reply())
    called = []
    turn = turns.take_turn(url=oauth.RESPONSES_URL,
                           headers={"Authorization": "Bearer synthetic-access"},
                           shape="responses", model="gpt-plan", system="",
                           messages=[{"role": "user", "text": "Which frame?"}], tools=[],
                           label="OpenAI", on_call=called.append)
    assert turn.text == "Frame 7 has the lowest RMSD."
    assert turn.calls == () and called == []
    assert (turn.usage.input_tokens, turn.usage.cache_read_tokens,
            turn.usage.output_tokens) == (10, 20, 6)
