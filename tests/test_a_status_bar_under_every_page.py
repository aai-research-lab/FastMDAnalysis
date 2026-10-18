"""A status bar under every page, as an editor's: a strip, not a panel.

It says this computer and its GPUs, the run going here with where it stands
and the time left, the open study's problems and the Agent's AI model, each
opening its page. The sidebar's progress card and two of the three Stop
buttons (the sidebar's, the builder's) went into it: one Stop, beside how
the run is going on its Overview.
"""

from __future__ import annotations

import json
import socket
from datetime import datetime, timezone
from pathlib import Path

import pytest

from tests.test_the_workspace_says_its_studies import _study

PAGES = ("overview", "viewer", "analysis", "report", "files", "studies", "run", "agent")

BAR = """() => {
    const bar = document.getElementById('status-bar').getBoundingClientRect();
    const box = (sel) => { const e = document.querySelector(sel); return e ? e.getBoundingClientRect() : null; };
    const main = box('.main'), sidebar = box('.sidebar');
    const shown = [...document.querySelectorAll('#status-bar .sb-item')]
        .filter((e) => e.offsetWidth > 0);
    return {
        top: Math.round(bar.top), bottom: Math.round(bar.bottom), height: Math.round(bar.height),
        width: Math.round(bar.width), mainBottom: Math.round(main.bottom),
        sidebarBottom: Math.round(sidebar.bottom), inner: innerHeight, innerWidth,
        items: Object.fromEntries(shown.map((e) => [e.id || e.className, e.textContent.trim()])),
        scrolls: document.getElementById('status-bar').scrollHeight
            > document.getElementById('status-bar').clientHeight + 1,
    };
}"""


def _serve(root: Path, look, *, width: int = 1440, hosted: bool = False):
    pytest.importorskip("playwright.sync_api")
    from playwright.sync_api import sync_playwright

    from fastmdxplora.gui.server import start_dashboard_session

    session = start_dashboard_session(output=str(root), host="127.0.0.1", port=0)
    try:
        with sync_playwright() as pw:
            browser = pw.chromium.launch(args=["--enable-unsafe-swiftshader"])
            page = browser.new_page(viewport={"width": width, "height": 900})
            page.set_default_timeout(60000)
            errors: list[str] = []
            page.on("pageerror", lambda error: errors.append(str(error)))
            said = look(page, session.url)
            browser.close()
    finally:
        session.server.shutdown()
    assert errors == []
    return said


def _running(root: Path) -> Path:
    (root / "simulation").mkdir(parents=True, exist_ok=True)
    (root / "simulation" / "live_status.json").write_text(json.dumps({
        "status": "running", "stage": "production", "current_step": 3000,
        "total_planned_steps": 12000, "timestep_fs": 2.0, "platform": "CPU",
        "simulation_time_completed_ns": 0.006,
        "last_update_timestamp": datetime.now(timezone.utc).isoformat(),
        "stage_states": {"setup": "completed", "minimization": "completed", "nvt": "completed",
                         "npt": "completed", "production": "running"}}), encoding="utf-8")
    return root


def test_the_strip_is_under_every_page_and_the_columns_end_above_it(tmp_path):
    study = _study(tmp_path / "ubiquitin", means={"rmsd": (0.1, 0.01)})

    def look(page, url):
        said = {}
        page.goto(url + "#overview", wait_until="domcontentloaded")
        page.wait_for_function("() => document.querySelector('#sb-machine:not([hidden])')")
        for name in PAGES:
            page.evaluate(f"() => window.FastMDXDashboard.navigate('{name}')")
            page.wait_for_timeout(600)
            said[name] = page.evaluate(BAR)
        page.set_viewport_size({"width": 390, "height": 844})
        page.wait_for_timeout(400)
        said["phone"] = page.evaluate(
            "() => getComputedStyle(document.getElementById('status-bar')).display")
        return said

    said = _serve(study, look)
    host = socket.gethostname().split(".")[0]
    for name in PAGES:
        bar = said[name]
        assert 22 <= bar["height"] <= 26 and bar["bottom"] == bar["inner"], name
        assert bar["width"] == bar["innerWidth"] and not bar["scrolls"], name
        assert bar["mainBottom"] <= bar["top"] and bar["sidebarBottom"] <= bar["top"], name
        assert bar["items"]["sb-machine"].startswith(host + " · "), name
        assert bar["items"]["sb-agent"] == "Agent · no AI model set", name
    assert said["phone"] == "none"


def test_a_run_going_is_in_the_strip_with_its_stage_percent_and_time_left(tmp_path):
    study = _running(_study(tmp_path / "going"))

    def look(page, url):
        page.goto(url + "#overview", wait_until="domcontentloaded")
        page.wait_for_function(
            "() => document.getElementById('status-run').dataset.run === 'running'")
        page.wait_for_timeout(500)
        run = page.eval_on_selector("#status-run", "r => [r.offsetWidth > 0, r.innerText, r.title]")
        stops = page.evaluate("""() => [...document.querySelectorAll('[data-stop-run]')]
            .map((b) => b.id)""")
        sidebar = page.evaluate("() => !!document.querySelector('.sidebar #sidebar-progress')")
        page.evaluate("() => window.FastMDXDashboard.navigate('analysis')")
        page.click("#status-run")
        page.wait_for_function("() => document.documentElement.dataset.page === 'overview'")
        return {"run": run, "stops": stops, "sidebar": sidebar}

    said = _serve(study, look)
    shown, text, title = said["run"]
    assert shown
    words = text.split()
    assert "Production" in words and "25.0%" in words and "left" in words
    assert "step 3000 of 12000" in title
    # One Stop, beside how the run is going.
    assert said["stops"] == ["health-stop"]
    assert said["sidebar"] is False


def test_the_machine_s_gpus_and_the_agent_s_model_are_said(tmp_path, monkeypatch):
    from fastmdxplora import gpu_here
    from fastmdxplora.agent.models import ModelChoice, save_choice
    from fastmdxplora.gui import status_bar
    from fastmdxplora.remote.gpu_room import Gpu, Room

    monkeypatch.setattr(gpu_here, "nvidia_smi_here", lambda: True)
    monkeypatch.setattr(gpu_here, "room_here", lambda: Room(gpus=(
        Gpu(index=0, uuid="GPU-1", name="NVIDIA GeForce RTX 4090", total_mb=24564,
            used_mb=1000, free_mb=23564, busy_pct=3),)))
    save_choice(ModelChoice(provider="anthropic", model="claude-test-model"))
    status_bar.forget()
    said = status_bar.status_bar()
    status_bar.forget()
    assert said["machine"]["gpus"][0]["name"] == "NVIDIA GeForce RTX 4090"
    assert said["machine"]["gpus_read"] is True
    assert said["agent"]["model"] == "claude-test-model"
    assert "key" not in json.dumps(said["agent"]).lower()

    study = _study(tmp_path / "ubiquitin")

    def look(page, url):
        page.goto(url + "#overview", wait_until="domcontentloaded")
        page.wait_for_function("() => document.querySelector('#sb-machine:not([hidden])')")
        page.wait_for_timeout(300)
        machine = page.eval_on_selector("#sb-machine", "b => [b.textContent.trim(), b.title]")
        agent = page.eval_on_selector("#sb-agent", "b => b.textContent.trim()")
        page.click("#sb-agent")
        page.wait_for_selector("#agent-settings:not([hidden])")
        return {"machine": machine, "agent": agent}

    status_bar.forget()
    seen = _serve(study, look)
    status_bar.forget()
    assert seen["machine"][0].endswith("· GeForce RTX 4090")
    assert "GPU 0: NVIDIA GeForce RTX 4090, 23 GB free of 24 GB, 3% busy" in seen["machine"][1]
    assert seen["agent"] == "Agent · claude-test-model"


def test_the_open_study_s_problems_are_counted_and_open_their_tab(tmp_path):
    study = _study(tmp_path / "failed", failed="simulation.run.unstable")

    def look(page, url):
        page.goto(url + "#overview", wait_until="domcontentloaded")
        page.wait_for_function("() => document.querySelector('#sb-problems:not([hidden])')")
        count = page.eval_on_selector("#sb-problems", "b => [b.dataset.count, b.textContent.trim()]")
        listed = page.evaluate("""async () => {
            const r = await fetch('/api/problems?study=' + encodeURIComponent(
                (await (await fetch('/api/status')).json()).study));
            return (await r.json()).items.length; }""")
        page.click("#sb-problems")
        page.wait_for_selector("#problems-pane:not([hidden])")
        return {"count": count, "listed": listed}

    said = _serve(study, look)
    assert said["count"][0] == str(said["listed"]) and said["listed"] >= 1
    assert said["count"][1].startswith(said["count"][0] + " problem")


def test_the_route_says_this_computer_here_and_nothing_to_a_service_s_guests(tmp_path):
    import threading
    import urllib.request
    from http.server import ThreadingHTTPServer

    from fastmdxplora.gui.exploration import DashboardRuntime
    from fastmdxplora.gui.hosting import SECRET_HEADER, Hosting
    from fastmdxplora.gui.server import GETS_ANSWERED_BEYOND_LOOPBACK, make_handler

    assert "/api/status-bar" not in GETS_ANSWERED_BEYOND_LOOPBACK
    study = _study(tmp_path / "person" / "study")
    secret = "s" * 40

    def ask(hosting):
        runtime = DashboardRuntime(workspace_root=study.parent, exploration_root=study.parent)
        runtime.hosting = hosting
        httpd = ThreadingHTTPServer(("127.0.0.1", 0), make_handler(
            study.parent, runtime=runtime, allow_control=True, hosting=hosting))
        threading.Thread(target=httpd.serve_forever, daemon=True).start()
        try:
            request = urllib.request.Request(
                f"http://127.0.0.1:{httpd.server_address[1]}/api/status-bar",
                headers={SECRET_HEADER: secret, "Host": "app.example.org"} if hosting else {})
            with urllib.request.urlopen(request, timeout=30) as response:
                return json.loads(response.read())
        finally:
            httpd.shutdown()
            httpd.server_close()

    here = ask(None)
    assert here["ok"] is True and here["machine"]["host"] == socket.gethostname().split(".")[0]
    assert here["machine"]["cpus"] >= 1
    hosted = ask(Hosting(workspace=study.parent.resolve(),
                         allowed_hosts=frozenset({"app.example.org"}), secret=secret))
    assert hosted["ok"] is False and "machine" not in hosted
