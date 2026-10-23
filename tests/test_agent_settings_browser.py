"""The Agent's provider settings are exercised in a real browser."""

from __future__ import annotations

import json

import pytest


def _stale_client_state(page) -> None:
    page.add_init_script("""
      localStorage.setItem("agent-auth", JSON.stringify({
        provider: "anthropic", auth_method: "anthropic_console", connected: true
      }));
      localStorage.setItem("agent-provider", "anthropic");
      sessionStorage.setItem("agent-account", JSON.stringify({
        provider: "openai", auth_method: "chatgpt_plan", signed_in: true
      }));
      history.replaceState({provider: "anthropic", auth_method: "anthropic_console",
                            connected: true}, "", location.href);
    """)


def _agent_page(page, url: str) -> None:
    page.goto(url + "#agent", wait_until="domcontentloaded")
    page.wait_for_function(
        "document.querySelectorAll('#agent-provider option').length === 2")


def test_agent_settings_offer_api_keys_without_auth_controls_or_restored_state(
        tmp_path, monkeypatch):
    pytest.importorskip("playwright.sync_api")
    from playwright.sync_api import sync_playwright

    from fastmdxplora.gui.server import start_dashboard_session

    settings = tmp_path / "settings"
    monkeypatch.setenv("FASTMDXPLORA_CONFIG_DIR", str(settings))
    monkeypatch.delenv("OPENAI_API_KEY", raising=False)
    monkeypatch.delenv("FASTMDX_MODEL_API_KEY", raising=False)
    session = start_dashboard_session(output=str(tmp_path), host="127.0.0.1", port=0)
    try:
        with sync_playwright() as playwright:
            browser = playwright.chromium.launch(args=["--enable-unsafe-swiftshader"])
            page = browser.new_page(viewport={"width": 1400, "height": 900})
            page.set_default_timeout(30000)
            _stale_client_state(page)
            _agent_page(page, session.url)

            assert page.locator("#agent-provider option").all_text_contents() == [
                "OpenAI", "Other (any OpenAI-compatible server)",
            ]
            assert page.locator("#agent-provider").input_value() == "openai"
            page.click("#agent-settings-open")
            page.wait_for_selector("#agent-settings:not([hidden])")
            assert page.locator("#agent-settings select").count() == 3
            for selector in (
                "#agent-auth-method", "#agent-auth-connect", "#agent-account",
                "#agent-account-add", "#agent-account-use", "#agent-auth-test",
                "#agent-auth-disconnect", "#agent-cache-for",
            ):
                assert page.locator(selector).count() == 0
            settings_text = page.locator("#agent-settings").inner_text().lower()
            for label in (
                "authentication", "continue with", "add account", "use account",
                "test connection", "disconnect",
            ):
                assert label not in settings_text
            assert "connected" not in settings_text and "signed in" not in settings_text
            assert "no ai model is set" in settings_text
            browser.close()
    finally:
        session.server.shutdown()


def test_saving_an_openai_key_is_shown_after_reload_without_echoing_it(
        tmp_path, monkeypatch):
    pytest.importorskip("playwright.sync_api")
    from playwright.sync_api import sync_playwright

    from fastmdxplora.agent import models
    from fastmdxplora.gui.server import start_dashboard_session

    settings = tmp_path / "settings"
    monkeypatch.setenv("FASTMDXPLORA_CONFIG_DIR", str(settings))
    monkeypatch.delenv("OPENAI_API_KEY", raising=False)
    monkeypatch.delenv("FASTMDX_MODEL_API_KEY", raising=False)
    monkeypatch.setattr(models, "list_model_options",
                        lambda _choice: [{"id": "gpt-5", "label": "gpt-5"}])
    session = start_dashboard_session(output=str(tmp_path), host="127.0.0.1", port=0)
    try:
        with sync_playwright() as playwright:
            browser = playwright.chromium.launch(args=["--enable-unsafe-swiftshader"])
            page = browser.new_page(viewport={"width": 1400, "height": 900})
            page.set_default_timeout(30000)
            _stale_client_state(page)
            model_responses = []

            def record_response(response):
                if "/api/agent/model" in response.url:
                    model_responses.append(response.text())

            page.on("response", record_response)
            _agent_page(page, session.url)
            page.click("#agent-settings-open")
            page.wait_for_selector("#agent-settings:not([hidden])")
            page.select_option("#agent-provider", "openai")
            page.select_option("#agent-model", "gpt-5")
            page.fill("#agent-key", "browser-test-key")
            page.click("#agent-save-model")
            page.wait_for_function(
                "document.querySelector('#agent-key-help').textContent.includes("
                "'An API key is set')")
            assert page.locator("#agent-key").input_value() == ""
            assert "browser-test-key" not in page.locator("body").inner_text()

            page.reload(wait_until="domcontentloaded")
            page.wait_for_function(
                "document.querySelector('#agent-key-help').textContent.includes("
                "'An API key is set')")
            assert page.locator("#agent-provider").input_value() == "openai"
            assert page.locator("#agent-key").input_value() == ""
            assert "browser-test-key" not in page.locator("body").inner_text()
            assert all("browser-test-key" not in body for body in model_responses)
            browser.close()
    finally:
        session.server.shutdown()


def test_a_stale_model_file_shows_no_model_or_signed_in_state_in_the_browser(
        tmp_path, monkeypatch):
    pytest.importorskip("playwright.sync_api")
    from playwright.sync_api import sync_playwright

    from fastmdxplora.gui.server import start_dashboard_session

    settings = tmp_path / "settings"
    settings.mkdir()
    (settings / "model.json").write_text(json.dumps({
        "provider": "anthropic",
        "model": "old-model",
        "cache_for": "one-hour",
        "auth_method": "chatgpt_plan",
        "oauth_access_token": "stale-browser-test-token",
        "account_id": "stale-account",
        "session": {"id": "stale-session"},
    }), encoding="utf-8")
    monkeypatch.setenv("FASTMDXPLORA_CONFIG_DIR", str(settings))
    monkeypatch.delenv("OPENAI_API_KEY", raising=False)
    monkeypatch.delenv("FASTMDX_MODEL_API_KEY", raising=False)
    session = start_dashboard_session(output=str(tmp_path), host="127.0.0.1", port=0)
    try:
        with sync_playwright() as playwright:
            browser = playwright.chromium.launch(args=["--enable-unsafe-swiftshader"])
            page = browser.new_page(viewport={"width": 1400, "height": 900})
            page.set_default_timeout(30000)
            _stale_client_state(page)
            model_responses = []
            page.on("response", lambda response: model_responses.append(response.text())
                    if "/api/agent/model" in response.url else None)
            _agent_page(page, session.url)
            page.click("#agent-settings-open")
            page.wait_for_selector("#agent-settings:not([hidden])")

            assert page.locator("#agent-provider").input_value() == "openai"
            assert "no ai model is set" in page.locator("#agent-model-summary").inner_text().lower()
            settings_text = page.locator("#agent-settings").inner_text().lower()
            assert "connected" not in settings_text and "signed in" not in settings_text
            assert "stale-browser-test-token" not in page.locator("body").inner_text()
            assert all("stale-browser-test-token" not in body for body in model_responses)
            browser.close()
    finally:
        session.server.shutdown()
