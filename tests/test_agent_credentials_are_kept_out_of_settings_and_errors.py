"""The provider route, its credential and its failures, outside the happy path.

``fastmdxplora.agent.models`` keeps API keys in the OS credential vault and
the settings file free of them, and turns each provider's failure into a
message Settings can show. ``fastmdxplora.agent.anthropic_cli`` drives the
official ``ant`` CLI in an app-owned profile. These tests cover the refusals:
an insecure or locked vault, a legacy plaintext key, an authentication method
a provider does not offer, every way a connection test can fail, and every
way the CLI can be missing, slow or fail. Every key and token is synthetic,
no request leaves the machine, and no real terminal or CLI is started.
"""

from __future__ import annotations

import importlib.util
import json
import subprocess
import sys
import types
import urllib.error
import urllib.request
from io import BytesIO
from subprocess import CompletedProcess

import pytest

from fastmdxplora.refusals import BackendUnavailable, StudyError


@pytest.fixture
def config(monkeypatch, tmp_path):
    monkeypatch.setenv("FASTMDXPLORA_CONFIG_DIR", str(tmp_path))
    for name in ("ANTHROPIC_API_KEY", "OPENAI_API_KEY", "FASTMDX_MODEL_API_KEY"):
        monkeypatch.delenv(name, raising=False)
    return tmp_path


class _Vault:
    """A synthetic keyring that can be told to refuse."""

    def __init__(self, *, refuse=()):
        self.values = {}
        self.refuse = set(refuse)

    def get_password(self, service, username):
        if "get" in self.refuse:
            raise RuntimeError("synthetic locked vault")
        return self.values.get((service, username))

    def set_password(self, service, username, value):
        if "set" in self.refuse:
            raise RuntimeError("synthetic denied write")
        self.values[(service, username)] = value

    def delete_password(self, service, username):
        if "delete" in self.refuse:
            raise RuntimeError("synthetic denied delete")
        del self.values[(service, username)]


def _use_vault(monkeypatch, vault):
    from fastmdxplora.agent import models

    monkeypatch.setattr(models, "_secure_keyring", lambda: vault)
    return vault


def _unpatched_secure_keyring():
    """The real vault selector; the suite's autouse fixture replaces it on the module."""
    spec = importlib.util.find_spec("fastmdxplora.agent.models")
    fresh = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(fresh)
    return fresh._secure_keyring


# -- the settings file and the vault ------------------------------------------


def test_console_and_plan_sign_in_send_their_token_as_a_bearer():
    from fastmdxplora.agent.models import ModelChoice

    assert ModelChoice("anthropic", "m").auth_style == "x-api-key"
    assert ModelChoice("anthropic", "m", auth_method="anthropic_console").auth_style == "bearer"
    assert ModelChoice("openai", "m", auth_method="chatgpt_plan").auth_style == "bearer"


def test_a_settings_file_that_is_not_an_object_means_nothing_is_chosen(config):
    from fastmdxplora.agent import models

    (config / "model.json").write_text(json.dumps(["openai", "gpt-5"]), encoding="utf-8")
    assert models._settings_record() == {}
    assert models.load_choice() is None


def test_a_legacy_plaintext_key_is_removed_and_never_returned(config):
    from fastmdxplora.agent import models

    path = config / "model.json"
    path.write_text(json.dumps({"provider": "openai", "model": "gpt-5",
                                "api_key": "sk-synthetic-legacy"}), encoding="utf-8")
    choice = models.load_choice()
    assert (choice.provider, choice.model, choice.auth_method) == ("openai", "gpt-5", "api_key")
    assert "sk-synthetic-legacy" not in path.read_text(encoding="utf-8")
    assert json.loads(path.read_text(encoding="utf-8")) == {"provider": "openai",
                                                            "model": "gpt-5"}


def test_an_authentication_method_the_provider_does_not_offer_falls_back_to_a_key(config):
    from fastmdxplora.agent import models

    (config / "model.json").write_text(json.dumps({
        "provider": "anthropic", "model": "m", "auth_method": "chatgpt_plan"}))
    assert models.load_choice().auth_method == "api_key"


@pytest.mark.parametrize("module,secure", [
    ("keyring.backends.Windows", True),
    ("keyring.backends.macOS.api", True),
    ("keyring.backends.SecretService", True),
    ("keyring.backends.libsecret", True),
    ("keyring.backends.fail", False),
    ("keyrings.alt.file", False),
])
def test_only_an_os_backed_keyring_is_trusted_with_keys(monkeypatch, module, secure):
    backend = type("Backend", (), {"__module__": module})()
    keyring = types.ModuleType("keyring")
    keyring.get_keyring = lambda: backend
    monkeypatch.setitem(sys.modules, "keyring", keyring)
    assert (_unpatched_secure_keyring()() is keyring) is secure


def test_an_unconfigured_keyring_is_no_vault(monkeypatch):
    keyring = types.ModuleType("keyring")

    def unconfigured():
        raise RuntimeError("synthetic: no recommended backend")

    keyring.get_keyring = unconfigured
    monkeypatch.setitem(sys.modules, "keyring", keyring)
    assert _unpatched_secure_keyring()() is None


def test_without_a_vault_keys_are_neither_read_saved_nor_forgotten(monkeypatch, config):
    from fastmdxplora.agent import models

    monkeypatch.setattr(models, "_secure_keyring", lambda: None)
    choice = models.ModelChoice("openai", "gpt-5")
    assert models._stored_api_key(choice) == ""
    with pytest.raises(StudyError) as caught:
        models.save_choice(choice, key="sk-synthetic")
    assert caught.value.code == "environment.credentials.secure_store_unavailable"
    assert not (config / "model.json").exists()
    with pytest.raises(StudyError, match="credential vault is unavailable") as caught:
        models.forget_api_key(choice)
    assert caught.value.code == "environment.credentials.secure_store_unavailable"


def test_a_locked_vault_reads_as_no_key_and_refuses_to_store_or_remove(monkeypatch, config):
    from fastmdxplora.agent import models

    choice = models.ModelChoice("openai", "gpt-5")
    _use_vault(monkeypatch, _Vault(refuse={"get", "set"}))
    assert models._stored_api_key(choice) == ""
    assert "nowhere" in models.describe_choice(models.save_choice(choice))
    with pytest.raises(StudyError, match="could not store this key") as caught:
        models.save_choice(choice, key="sk-synthetic-denied")
    assert "sk-synthetic-denied" not in str(caught.value)
    with pytest.raises(StudyError, match="could not be removed"):
        models.forget_api_key(choice)


def test_forgetting_a_key_that_was_never_saved_reports_nothing_removed(monkeypatch, config):
    from fastmdxplora.agent import models

    vault = _use_vault(monkeypatch, _Vault())
    choice = models.ModelChoice("openai", "gpt-5")
    assert models.forget_api_key(choice) is False
    models.save_choice(choice, key="sk-synthetic")
    assert models.forget_api_key(choice) is True and vault.values == {}


def test_a_route_is_saved_only_with_an_offered_method_and_a_key_only_for_api_keys(config):
    from fastmdxplora.agent import models

    with pytest.raises(StudyError, match="not supported for this provider") as caught:
        models.save_choice(models.ModelChoice("compatible", "m", auth_method="chatgpt_plan"))
    assert caught.value.code == "config.option.not_permitted"
    with pytest.raises(StudyError, match="cannot be saved for this authentication") as caught:
        models.save_choice(models.ModelChoice("openai", "gpt-5", auth_method="chatgpt_plan"),
                           key="sk-synthetic")
    assert caught.value.code == "config.option.not_permitted"
    assert not (config / "model.json").exists()


@pytest.mark.parametrize("method,module,error", [
    ("anthropic_console", "anthropic_cli", "ConsoleSignInError"),
    ("anthropic_console", "anthropic_cli", "BackendUnavailable"),
    ("chatgpt_plan", "openai_chatgpt", "SignInError"),
])
def test_a_sign_in_failure_keeps_its_code_when_a_key_is_asked_for(
        monkeypatch, config, method, module, error):
    import importlib

    from fastmdxplora.agent import models

    provider = "anthropic" if module == "anthropic_cli" else "openai"
    source = importlib.import_module(f"fastmdxplora.agent.{module}")
    kind = BackendUnavailable if error == "BackendUnavailable" else getattr(source, error)
    failure = kind("Synthetic sign-in failure.")

    def refuse():
        raise failure

    monkeypatch.setattr(source, "access_token", refuse)
    with pytest.raises(StudyError, match="Synthetic sign-in failure") as caught:
        models._key_for(models.ModelChoice(provider, "m", auth_method=method))
    assert caught.value.code == failure.code
    assert caught.value.refusal.details["provider"] == provider


def test_each_sign_in_route_is_described_by_how_it_authenticates(config):
    from fastmdxplora.agent import models

    models.save_choice(models.ModelChoice("anthropic", "m", auth_method="anthropic_console"))
    assert "Anthropic Console OAuth through the official ant CLI" in models.describe_choice()
    models.save_choice(models.ModelChoice("openai", "gpt-5", auth_method="chatgpt_plan"))
    said = models.describe_choice()
    assert "Sign in with ChatGPT plan usage" in said and "API key" not in said


# -- listing and testing the route ------------------------------------------


class _Reply(BytesIO):
    def __enter__(self):
        return self

    def __exit__(self, *_args):
        self.close()


def _serve(monkeypatch, reply):
    asked = []

    def urlopen(request, timeout=None):
        asked.append(request)
        if isinstance(reply, BaseException):
            raise reply
        return _Reply(reply)

    monkeypatch.setattr(urllib.request, "urlopen", urlopen)
    return asked


def test_plan_models_come_from_the_plan_catalog_or_the_fallback(monkeypatch, config):
    from fastmdxplora.agent import models
    from fastmdxplora.agent import openai_chatgpt as oauth

    choice = models.ModelChoice("openai", "gpt-5", auth_method="chatgpt_plan")
    monkeypatch.setattr(oauth, "list_models", lambda timeout: [
        {"id": "gpt-plan-a", "label": "A"}, {"id": "gpt-plan-b", "label": "B"}])
    assert models.list_models(choice) == ("gpt-plan-a", "gpt-plan-b")
    assert models.list_model_options(choice)[1] == {"id": "gpt-plan-b", "label": "B"}
    monkeypatch.setattr(oauth, "list_models", lambda timeout: (_ for _ in ()).throw(
        oauth.SignInError("No active ChatGPT sign-in.")))
    assert models.list_models(choice) == models.PROVIDERS["openai"]["models"]


def test_a_route_is_tested_only_once_it_is_saved(config):
    from fastmdxplora.agent import models

    assert models.probe_connection() == {"connected": False,
                                         "message": "Save a provider and model first."}


def test_a_plan_route_test_reports_sign_in_and_model_availability(monkeypatch, config):
    from fastmdxplora.agent import models
    from fastmdxplora.agent import openai_chatgpt as oauth

    options = [{"id": "gpt-plan", "label": "GPT plan"}]
    monkeypatch.setattr(oauth, "list_models", lambda timeout: options)
    verified = models.probe_connection(
        models.ModelChoice("openai", "gpt-plan", auth_method="chatgpt_plan"))
    assert verified["connected"] is True
    assert verified["message"] == "ChatGPT plan route verified. 1 available model(s)."
    assert verified["model_options"] == options
    missing = models.probe_connection(
        models.ModelChoice("openai", "gpt-unlisted", auth_method="chatgpt_plan"))
    assert missing["connected"] is False and missing["models"] == ["gpt-plan"]
    assert "not available to this account" in missing["message"]
    monkeypatch.setattr(oauth, "list_models", lambda timeout: (_ for _ in ()).throw(
        oauth.SignInError("No active ChatGPT sign-in.")))
    refused = models.probe_connection(
        models.ModelChoice("openai", "gpt-plan", auth_method="chatgpt_plan"))
    assert refused == {"connected": False, "message": "No active ChatGPT sign-in.",
                       "code": "environment.credentials.sign_in_failed"}


def test_a_compatible_server_without_a_base_url_has_nothing_to_test(config):
    from fastmdxplora.agent import models

    answer = models.probe_connection(models.ModelChoice("compatible", "local-model"))
    assert answer == {"connected": False,
                      "message": "This provider has no model-list endpoint to test."}


def test_a_route_without_a_key_is_reported_before_any_request(monkeypatch, config):
    from fastmdxplora.agent import models

    asked = _serve(monkeypatch, b"{}")
    answer = models.probe_connection(models.ModelChoice("openai", "gpt-5"))
    assert answer["connected"] is False
    assert answer["code"] == "environment.credentials.absent"
    assert "OPENAI_API_KEY" in answer["message"] and asked == []


@pytest.mark.parametrize("payload,connected,message", [
    ({"data": [{"id": "gpt-5"}, {"id": "gpt-other"}, {"name": "no-id"}]}, True,
     "Connected. 2 available model(s) reported by the provider."),
    ({"data": [{"id": "gpt-other"}]}, True, "selected model is not visible"),
    ({"data": "gpt-5"}, False, "The provider response was not a model list."),
    (["gpt-5"], False, "The provider response was not a model list."),
])
def test_a_key_route_test_reads_the_provider_model_list(monkeypatch, config, payload,
                                                         connected, message):
    from fastmdxplora.agent import models

    monkeypatch.setenv("FASTMDX_MODEL_API_KEY", "sk-synthetic-compatible")
    asked = _serve(monkeypatch, json.dumps(payload).encode())
    answer = models.probe_connection(models.ModelChoice(
        "compatible", "gpt-5", base_url="https://compatible.example.invalid/v1/"))
    assert answer["connected"] is connected and message in answer["message"]
    assert asked[0].full_url == "https://compatible.example.invalid/v1/models"
    assert asked[0].get_header("Authorization") == "Bearer sk-synthetic-compatible"
    assert "sk-synthetic" not in json.dumps(answer)


@pytest.mark.parametrize("method,failure,message", [
    ("api_key", 401, "The API key was rejected."),
    ("anthropic_console", 401, "Anthropic Console login expired"),
    ("api_key", 403, "does not have access to this provider resource"),
    ("api_key", 429, "rate limit or account quota"),
    ("api_key", 500, "returned HTTP 500"),
    ("api_key", urllib.error.URLError("synthetic"), "could not be reached"),
    ("api_key", TimeoutError(), "could not be reached"),
    ("api_key", b"not json", "unreadable response"),
])
def test_a_failed_route_test_says_what_failed_without_the_credential(
        monkeypatch, config, method, failure, message):
    from fastmdxplora.agent import anthropic_cli, models

    monkeypatch.setenv("ANTHROPIC_API_KEY", "sk-ant-synthetic")
    monkeypatch.setattr(anthropic_cli, "access_token", lambda: "console-synthetic")
    if isinstance(failure, int):
        failure = urllib.error.HTTPError(models.PROVIDERS["anthropic"]["models_url"],
                                         failure, "synthetic", {}, BytesIO(b"sk-ant-synthetic"))
    asked = _serve(monkeypatch, failure)
    answer = models.probe_connection(models.ModelChoice("anthropic", "m", auth_method=method))
    assert answer["connected"] is False and message in answer["message"]
    assert "synthetic" not in answer["message"]
    expected = {"x-api-key": "sk-ant-synthetic"} if method == "api_key" else {
        "Authorization": "Bearer console-synthetic"}
    for name, value in expected.items():
        assert asked[0].get_header(name.capitalize()) == value


def test_a_plan_completion_is_one_responses_turn_with_the_plan_token(monkeypatch, config):
    from fastmdxplora.agent import models, turns
    from fastmdxplora.agent import openai_chatgpt as oauth

    monkeypatch.setattr(oauth, "access_token", lambda: "synthetic-plan-token")
    taken = []

    def take_turn(**kwargs):
        taken.append(kwargs)
        kwargs["on_text"]("streamed ")
        return types.SimpleNamespace(text="streamed reply")

    monkeypatch.setattr(turns, "take_turn", take_turn)
    pieces = []
    complete = models.completion_for(
        models.ModelChoice("openai", "gpt-plan", auth_method="chatgpt_plan"), timeout=9)
    assert complete("What is RMSD?", pieces.append) == "streamed reply"
    assert pieces == ["streamed "]
    turn = taken[0]
    assert turn["url"] == oauth.RESPONSES_URL and turn["shape"] == "responses"
    assert turn["headers"]["Authorization"] == "Bearer synthetic-plan-token"
    assert turn["messages"] == [{"role": "user", "text": "What is RMSD?"}]
    assert (turn["model"], turn["tools"], turn["timeout"]) == ("gpt-plan", [], 9)


# -- the Anthropic CLI profile ------------------------------------------------


@pytest.fixture
def cli(monkeypatch, tmp_path):
    from fastmdxplora.agent import anthropic_cli

    monkeypatch.setattr(anthropic_cli, "user_config_dir", lambda: tmp_path)
    monkeypatch.setattr(anthropic_cli, "executable", lambda: "/opt/synthetic/ant")
    monkeypatch.setattr(anthropic_cli, "_login_pending", False)
    monkeypatch.setattr(anthropic_cli, "_login_process", None)
    monkeypatch.setattr(anthropic_cli, "_login_started", None)
    return anthropic_cli


def _platform(monkeypatch, cli, *, os_name="posix", platform="linux", which=None):
    monkeypatch.setattr(cli, "os", types.SimpleNamespace(
        name=os_name, environ={"PATH": "/usr/bin", "HOME": "/home/synthetic"}))
    monkeypatch.setattr(cli, "sys", types.SimpleNamespace(platform=platform))
    found = dict(which or {})
    monkeypatch.setattr(cli.shutil, "which", lambda name: found.get(name))


@pytest.mark.parametrize("profile,message", [
    ({"generation": "not-hex"}, "identifier is invalid"),
    ({"generation": "A" * 32}, "identifier is invalid"),
    ({"generation": 7}, "identifier is invalid"),
    ({"other": "0" * 32}, "could not be read"),
    ("{broken", "could not be read"),
])
def test_an_unreadable_console_profile_record_is_refused(cli, tmp_path, profile, message):
    raw = profile if isinstance(profile, str) else json.dumps(profile)
    (tmp_path / "anthropic-console-profile.json").write_text(raw, encoding="utf-8")
    with pytest.raises(cli.ConsoleSignInError, match=message):
        cli.config_dir()


def test_the_console_profile_is_the_first_one_until_it_is_retired(cli, tmp_path):
    assert cli.config_dir() == tmp_path / "anthropic-console"
    cli._retire_profile()
    generation = json.loads((tmp_path / "anthropic-console-profile.json").read_text())[
        "generation"]
    assert cli.config_dir() == tmp_path / f"anthropic-console-{generation}"


def test_a_disconnect_that_cannot_retire_the_profile_is_refused(monkeypatch, cli):
    from fastmdxplora.agent import models

    ran = []
    monkeypatch.setattr(cli, "_run", lambda *a, **kw: ran.append(a))

    def read_only(*_args):
        raise PermissionError("synthetic read-only settings")

    monkeypatch.setattr(models, "_write_settings", read_only)
    with pytest.raises(cli.ConsoleSignInError, match="disconnect could not be saved"):
        cli.disconnect()
    assert ran == []


@pytest.mark.parametrize("os_name,platform,which,expected", [
    ("posix", "linux", {}, False),
    ("nt", "win32", {"ant": "C:/ant.exe"}, True),
    ("posix", "darwin", {"ant": "/ant"}, False),
    ("posix", "darwin", {"ant": "/ant", "osascript": "/usr/bin/osascript"}, True),
    ("posix", "linux", {"ant": "/ant"}, False),
    ("posix", "linux", {"ant": "/ant", "gnome-terminal": "/usr/bin/gnome-terminal"}, True),
    ("posix", "linux", {"ant": "/ant", "x-terminal-emulator": "/usr/bin/xte"}, True),
])
def test_login_is_offered_only_where_a_visible_terminal_can_open(
        monkeypatch, cli, os_name, platform, which, expected):
    from fastmdxplora.agent import anthropic_cli

    _platform(monkeypatch, cli, os_name=os_name, platform=platform, which=which)
    monkeypatch.setattr(cli, "executable", lambda: anthropic_cli.shutil.which("ant"))
    assert cli.can_launch_login() is expected


def test_without_the_cli_nothing_is_run_and_each_entry_point_says_so(monkeypatch, cli):
    monkeypatch.setattr(cli, "executable", lambda: None)
    monkeypatch.setattr(cli.subprocess, "run", lambda *a, **kw: pytest.fail("ran ant"))
    with pytest.raises(BackendUnavailable, match="not installed"):
        cli._run(["models", "list"])
    with pytest.raises(BackendUnavailable, match="Install Anthropic's official ant CLI"):
        cli.access_token()
    assert cli.status() == {
        "available": False, "connected": False,
        "message": "Install Anthropic's official ant CLI to use Anthropic Console sign-in."}
    started = cli.connect()
    assert started["ok"] is False and "Install Anthropic's official ant CLI" in started["message"]
    assert cli.disconnect() == {
        "ok": True, "connected": False,
        "message": "This app has no Anthropic Console profile to disconnect."}


def test_the_cli_runs_with_the_app_profile_and_no_terminal_input(monkeypatch, cli, tmp_path):
    seen = []
    monkeypatch.setattr(cli.subprocess, "run", lambda args, **kw: seen.append((args, kw)) or
                        CompletedProcess(args, 0, stdout="", stderr=""))
    cli._run(["models", "list"], timeout=4, profile=tmp_path / "chosen")
    args, kwargs = seen[0]
    assert args == ["/opt/synthetic/ant", "models", "list"]
    assert kwargs["env"]["ANTHROPIC_CONFIG_DIR"] == str(tmp_path / "chosen")
    assert kwargs["stdin"] is subprocess.DEVNULL and kwargs["timeout"] == 4
    assert kwargs["check"] is False and kwargs["capture_output"] is True


def test_status_does_not_confirm_a_profile_retired_while_it_was_checked(monkeypatch, cli):
    def retire_during_check(args, **_kwargs):
        cli._retire_profile()
        return CompletedProcess(args, 0, stdout="[]", stderr="")

    monkeypatch.setattr(cli, "_run", retire_during_check)
    assert cli.status() == {"available": True, "connected": False,
                            "message": "Anthropic Console sign-in changed. Refresh its status."}


@pytest.mark.parametrize("failure,message", [
    (subprocess.TimeoutExpired("ant", 15), "did not respond"),
    (OSError("synthetic"), "could not be started"),
])
def test_a_slow_or_unstartable_cli_is_not_connected(monkeypatch, cli, failure, message):
    def fail(*_args, **_kwargs):
        raise failure

    monkeypatch.setattr(cli, "_run", fail)
    answer = cli.status()
    assert answer["connected"] is False and message in answer["message"]
    with pytest.raises(cli.ConsoleSignInError) as caught:
        cli.access_token()
    assert ("timed out" if isinstance(failure, subprocess.TimeoutExpired)
            else "could not be started") in str(caught.value)


def test_an_abandoned_login_stops_waiting_after_five_minutes(monkeypatch, cli):
    monkeypatch.setattr(cli, "_run", lambda args, **_kw: CompletedProcess(args, 1, "", ""))
    monkeypatch.setattr(cli, "_login_pending", True)
    monkeypatch.setattr(cli, "_login_process", None)
    monkeypatch.setattr(cli, "_login_started", cli.time.monotonic() - 301)
    answer = cli.status()
    assert answer["connected"] is False and "canceled, failed, or timed out" in answer["message"]
    assert cli._login_pending is False


def test_an_empty_token_is_not_a_sign_in(monkeypatch, cli):
    monkeypatch.setattr(cli, "_run", lambda args, **_kw: CompletedProcess(args, 0, "  \n", ""))
    with pytest.raises(cli.ConsoleSignInError, match="expired or is unavailable"):
        cli.access_token()


def test_only_one_login_waits_at_a_time(monkeypatch, cli):
    monkeypatch.setattr(cli, "_login_pending", True)
    monkeypatch.setattr(cli.subprocess, "Popen", lambda *a, **kw: pytest.fail("opened twice"))
    answer = cli.connect()
    assert answer["ok"] is False and "already waiting" in answer["message"]


def test_windows_login_opens_in_a_new_console_with_the_app_profile(monkeypatch, cli, tmp_path):
    _platform(monkeypatch, cli, os_name="nt", platform="win32")
    monkeypatch.setattr(cli.subprocess, "CREATE_NEW_CONSOLE", 0x10, raising=False)
    launched = []
    monkeypatch.setattr(cli.subprocess, "Popen",
                        lambda args, **kw: launched.append((args, kw)) or "process")
    assert cli.connect()["ok"] is True
    args, kwargs = launched[0]
    assert args == ["/opt/synthetic/ant", "auth", "login"]
    assert kwargs["creationflags"] == 0x10
    assert kwargs["env"]["ANTHROPIC_CONFIG_DIR"] == str(tmp_path / "anthropic-console")
    assert cli._login_pending is True and cli._login_process == "process"


@pytest.mark.parametrize("which,expected", [
    ({"x-terminal-emulator": "/usr/bin/x-terminal-emulator"},
     ["/usr/bin/x-terminal-emulator", "-e", "/opt/synthetic/ant", "auth", "login"]),
    ({"gnome-terminal": "/usr/bin/gnome-terminal"},
     ["/usr/bin/gnome-terminal", "--wait", "--", "/opt/synthetic/ant", "auth", "login"]),
    ({}, ["gnome-terminal", "--wait", "--", "/opt/synthetic/ant", "auth", "login"]),
])
def test_linux_login_opens_in_the_system_terminal(monkeypatch, cli, which, expected):
    _platform(monkeypatch, cli, which=which)
    launched = []
    monkeypatch.setattr(cli.subprocess, "Popen",
                        lambda args, **kw: launched.append((args, kw)) or "process")
    assert cli.connect()["ok"] is True
    args, kwargs = launched[0]
    assert args == expected
    assert kwargs["stdin"] is subprocess.DEVNULL and kwargs["close_fds"] is True
    assert "ANTHROPIC_API_KEY" not in kwargs["env"]


def test_a_terminal_that_cannot_open_leaves_no_login_waiting(monkeypatch, cli):
    _platform(monkeypatch, cli)

    def no_terminal(*_args, **_kwargs):
        raise FileNotFoundError("gnome-terminal")

    monkeypatch.setattr(cli.subprocess, "Popen", no_terminal)
    answer = cli.connect()
    assert answer["ok"] is False and "terminal could not be opened" in answer["message"]
    assert cli._login_pending is False


class _Login:
    def __init__(self, *, stubborn=False):
        self.stubborn = stubborn
        self.calls = []

    def poll(self):
        return None

    def terminate(self):
        self.calls.append("terminate")

    def wait(self, timeout):
        self.calls.append(("wait", timeout))
        if self.stubborn:
            raise subprocess.TimeoutExpired("ant", timeout)

    def kill(self):
        self.calls.append("kill")


@pytest.mark.parametrize("stubborn,calls", [
    (False, ["terminate", ("wait", 3)]),
    (True, ["terminate", ("wait", 3), "kill"]),
])
def test_disconnect_ends_a_waiting_login_terminal(monkeypatch, cli, stubborn, calls):
    login = _Login(stubborn=stubborn)
    monkeypatch.setattr(cli, "_login_pending", True)
    monkeypatch.setattr(cli, "_login_process", login)
    profiles = []
    monkeypatch.setattr(cli, "_run", lambda args, *, profile: profiles.append(profile) or
                        CompletedProcess(args, 0, "", ""))
    first = cli.config_dir()
    assert cli.disconnect()["message"] == "Disconnected from Anthropic Console."
    assert login.calls == calls
    assert cli._login_pending is False and cli._login_process is None
    # The logout runs against the profile the login was using, not the new one.
    assert profiles == [first] and cli.config_dir() != first


@pytest.mark.parametrize("outcome,message", [
    (subprocess.TimeoutExpired("ant", 15), "sign-out timed out"),
    (OSError("synthetic"), "could not be started"),
    (CompletedProcess(["ant"], 1, "", "synthetic-secret"), "sign-out did not complete"),
])
def test_a_failed_logout_is_not_reported_as_disconnected(monkeypatch, cli, outcome, message):
    def logout(args, **_kwargs):
        if isinstance(outcome, BaseException):
            raise outcome
        return outcome

    monkeypatch.setattr(cli, "_run", logout)
    answer = cli.disconnect()
    assert answer["ok"] is False and message in answer["message"]
    assert "connected" not in answer and "synthetic-secret" not in answer["message"]
