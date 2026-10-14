"""The GUI reaches your machines under the rules every interface keeps.

Its routes (`gui/remote_routes.py`) do what `fastmdx remote` does at a
terminal, for the page: machines only from their records, a send only of
the plan the page showed (a token kept ten minutes, used once, and only
while what would travel is what was shown), a fetch only of the size shown,
everything inside the GUI's workspace, behind the server's gate, and none
of it in a hosted GUI. The machine is this computer standing in, as in
`test_a_study_travels_and_comes_back`.
"""

from __future__ import annotations

import contextlib
import json
import threading
import time
import urllib.error
import urllib.request
from http.server import ThreadingHTTPServer
from pathlib import Path

import pytest

from fastmdxplora.gui.exploration import DashboardRuntime
from fastmdxplora.gui.remote_routes import RemoteDesk
from fastmdxplora.gui.server import make_handler
from tests import test_a_study_travels_and_comes_back as travels
from tests.test_a_study_travels_and_comes_back import RELEASE

machine = travels.machine
machine_path = travels.machine_path


@pytest.fixture
def served(machine, monkeypatch):
    """A GUI on the study's folder, reaching the stand-in machine."""
    from fastmdxplora.remote import api
    from fastmdxplora.remote import send as sending

    monkeypatch.setattr(api, "_link", lambda name, transport: machine.transport())
    monkeypatch.setattr(api, "this_code", lambda: RELEASE)
    monkeypatch.setattr(sending, "this_code", lambda: RELEASE)
    real = sending.run_here
    monkeypatch.setattr(sending, "run_here", lambda command, runner=None, **more:
                        real(command, runner=machine.local, **more))
    root = machine.study.parent
    # Started in the folder it opens, as `fastmdx gui` is: a relative config
    # name is the same file to the check and to the plan.
    monkeypatch.chdir(root)
    with _serving(DashboardRuntime(workspace_root=root, exploration_root=root)) as address:
        yield address, machine


@contextlib.contextmanager
def _serving(runtime: DashboardRuntime, *, allow_control: bool = True):
    handler = make_handler(runtime.workspace_root, runtime=runtime,
                           allow_control=allow_control)
    httpd = ThreadingHTTPServer(("127.0.0.1", 0), handler)
    threading.Thread(target=httpd.serve_forever, daemon=True).start()
    try:
        yield f"127.0.0.1:{httpd.server_address[1]}"
    finally:
        httpd.shutdown()
        httpd.server_close()


def _ask(address: str, path: str, body: dict | None = None,
         headers: dict | None = None) -> tuple[int, dict]:
    data = None if body is None else json.dumps(body).encode()
    request = urllib.request.Request(
        f"http://{address}{path}", data=data, method="GET" if body is None else "POST",
        headers={"Content-Type": "application/json", **(headers or {})})
    try:
        with urllib.request.urlopen(request, timeout=60) as response:
            return response.getcode(), json.loads(response.read() or b"{}")
    except urllib.error.HTTPError as exc:
        return exc.code, json.loads(exc.read() or b"{}")


def _planned(address: str) -> dict:
    status, said = _ask(address, "/api/remote/plan", {"config": "study.yml",
                                                      "machine": "box"})
    assert status == 200 and said["ok"], said
    return said


def _sent(machine) -> bool:
    return any(c.startswith("mkdir") for c in machine.commands)


class TestMachinesAndPlans:
    def test_the_machines_are_listed_from_their_records(self, served):
        address, machine = served
        status, said = _ask(address, "/api/remote/machines")
        assert status == 200 and said["ok"]
        assert [m["name"] for m in said["machines"]] == ["box"]
        assert said["jobs"] == []
        assert machine.commands == []

    def test_a_plan_says_what_travels_and_sends_nothing(self, served):
        address, machine = served
        plan = _planned(address)
        assert plan["machine"] == "box" and plan["job"] == "study"
        assert plan["travels"] == [{"name": "inputs/top.pdb", "from": "top.pdb",
                                    "bytes": 5}]
        assert plan["results"] == "study" and plan["kept_s"] == 600
        assert "explore -c study.yml --output run" in plan["script"]
        assert not _sent(machine)

    def test_the_plan_shown_is_sent_once(self, served):
        address, machine = served
        plan = _planned(address)
        status, sent = _ask(address, "/api/remote/send", {"plan": plan["plan"]})
        assert status == 200 and sent["ok"], sent
        assert sent["job"]["name"] == "study" and sent["job"]["machine"] == "box"
        status, again = _ask(address, "/api/remote/send", {"plan": plan["plan"]})
        assert again["ok"] is False and again["code"] == "remote.send.unconfirmed"
        travels._until_finished(machine, "study")

    def test_a_config_changed_after_the_plan_is_not_sent(self, served):
        """The person agreed to what the page showed: a config changed since
        sends something else."""
        address, machine = served
        plan = _planned(address)
        machine.study.write_text(machine.study.read_text() + "simulation:\n"
                                 "  random_seed: 7\n")
        status, said = _ask(address, "/api/remote/send", {"plan": plan["plan"]})
        assert said["ok"] is False and said["code"] == "remote.send.unconfirmed"
        assert "changed after the plan was shown" in said["error"]
        assert not _sent(machine)

    def test_a_plan_older_than_ten_minutes_is_not_sent(self, served, monkeypatch):
        from fastmdxplora.gui import remote_routes

        address, machine = served
        plan = _planned(address)
        later = time.monotonic() + remote_routes.PLAN_KEPT_S + 1
        monkeypatch.setattr(remote_routes.time, "monotonic", lambda: later)
        status, said = _ask(address, "/api/remote/send", {"plan": plan["plan"]})
        assert said["ok"] is False and said["code"] == "remote.send.unconfirmed"
        assert not _sent(machine)

    def test_a_token_not_given_out_sends_nothing(self, served):
        address, machine = served
        status, said = _ask(address, "/api/remote/send", {"plan": "0" * 32})
        assert said["ok"] is False and said["code"] == "remote.send.unconfirmed"
        assert machine.commands == []

    def test_a_config_outside_the_workspace_is_refused(self, served, tmp_path):
        address, machine = served
        root = machine.study.parent.resolve()
        far = (tmp_path / "far").resolve()
        far.mkdir()
        outside = far / "elsewhere.yml"
        outside.write_text(machine.study.read_text())
        (far / "top.pdb").write_text((root / "top.pdb").read_text())
        (far / "notes.txt").write_text("not a config\n")
        status, said = _ask(address, "/api/remote/plan", {"config": str(outside),
                                                          "machine": "box"})
        # Round 3 of the last round: it was called "not a config file", and
        # named no way to send it.
        assert said["ok"] is False and said["code"] == "remote.input.outside"
        assert said["error"].startswith(f"{outside} is outside the workspace")
        # Round 2: the copy named a file that may already be there.
        assert (f"Open the GUI on the config's folder (`cd {far}`, then `fastmdx gui` with "
                "no `--output`), or copy the config, with the files it names, into a new "
                f"folder inside {root} and give its full path there.") in said["error"]
        # Round 1: what is no config file is said as that, wherever it is;
        # remote-b4's round 3: with a way said alike for every such name.
        for named in (tmp_path / "gone.yml", far / "notes.txt"):
            _, none = _ask(address, "/api/remote/plan", {"config": str(named),
                                                         "machine": "box"})
            assert none["error"] == (
                f"{named} is not a config file in the workspace ({root}). To send a config "
                "kept elsewhere, copy it, with the files it names, into a new folder inside "
                f"{root} and give its full path there."), named
        # Each way, followed, plans it.
        (root / "new").mkdir()
        for name in ("elsewhere.yml", "top.pdb"):
            (root / "new" / name).write_bytes((far / name).read_bytes())
        _, copied = _ask(address, "/api/remote/plan", {
            "config": str(root / "new" / "elsewhere.yml"), "machine": "box"})
        assert copied["ok"], copied
        with _serving(DashboardRuntime(workspace_root=far, exploration_root=far.parent)
                      ) as there:
            _, opened = _ask(there, "/api/remote/plan", {"config": str(outside),
                                                         "machine": "box"})
        assert opened["ok"], opened
        assert machine.commands == [] or not _sent(machine)

    def test_a_results_folder_in_use_is_refused(self, served):
        address, machine = served
        (machine.study.parent / "study").mkdir()
        (machine.study.parent / "study" / "kept.txt").write_text("x")
        status, said = _ask(address, "/api/remote/plan", {"config": "study.yml",
                                                          "machine": "box"})
        assert said["ok"] is False and said["code"] == "environment.path.exists"

    def test_a_machine_without_a_record_is_refused(self, served):
        address, machine = served
        status, said = _ask(address, "/api/remote/plan", {"config": "study.yml",
                                                          "machine": "gpu-box"})
        assert said["ok"] is False and "gpu-box" in said["error"]
        assert machine.commands == []


class TestJobs:
    def _ended(self, served):
        address, machine = served
        plan = _planned(address)
        _ask(address, "/api/remote/send", {"plan": plan["plan"]})
        travels._until_finished(machine, "study")
        return address, machine

    def test_a_job_is_asked_about_and_fetched_at_the_size_shown(self, served):
        address, machine = self._ended(served)
        status, said = _ask(address, "/api/remote/job?job=study")
        assert said["ok"] and said["job"]["state"] == "done"
        status, sizes = _ask(address, "/api/remote/fetch-sizes?job=study")
        assert sizes["ok"] and sizes["run_written"]
        shown = sizes["bringing"]["without"]
        status, wrong = _ask(address, "/api/remote/fetch", {"job": "study",
                                                            "bringing": shown + 1})
        assert wrong["ok"] is False and wrong["code"] == "remote.fetch.unconfirmed"
        assert not (machine.study.parent / "study" / "manifest.json").exists()
        status, fetched = _ask(address, "/api/remote/fetch", {"job": "study",
                                                              "bringing": shown})
        assert fetched["ok"], fetched
        assert (machine.study.parent / "study" / "manifest.json").is_file()

    def test_a_job_of_another_workspace_is_not_reached(self, served):
        address, machine = served
        status, said = _ask(address, "/api/remote/job?job=trial")
        assert said["ok"] is False and said["code"] == "remote.job.unknown"

    def test_a_running_job_is_cancelled(self, served):
        address, machine = served
        machine.env["FAKE_SLEEP"] = "30"
        plan = _planned(address)
        status, sent = _ask(address, "/api/remote/send", {"plan": plan["plan"]})
        try:
            status, said = _ask(address, "/api/remote/cancel", {"job": "study"})
            assert said["ok"] and said["stopped"] and said["job"]["state"] == "abandoned"
        finally:
            from tests.test_remote_from_every_interface import _group_killed

            from fastmdxplora.remote.jobs import load_job
            _group_killed(load_job("study").handle)


class TestTheGate:
    @pytest.mark.parametrize("headers", [
        {"Origin": "https://evil.example"},
        {"Sec-Fetch-Site": "cross-site"},
    ])
    def test_a_page_elsewhere_sends_nothing(self, served, headers):
        address, machine = served
        plan = _planned(address)
        status, said = _ask(address, "/api/remote/send", {"plan": plan["plan"]},
                            headers=headers)
        assert status == 403
        assert not _sent(machine)

    def test_nothing_is_answered_beyond_loopback(self, machine):
        root = machine.study.parent
        with _serving(DashboardRuntime(workspace_root=root, exploration_root=root),
                      allow_control=False) as address:
            status, _ = _ask(address, "/api/remote/machines")
            assert status == 403
            status, _ = _ask(address, "/api/remote/plan", {"config": "study.yml",
                                                           "machine": "box"})
            assert status == 403
        assert machine.commands == []

    def test_a_hosted_gui_reaches_no_machine(self, machine):
        desk = RemoteDesk(machine.study.parent, hosted=True)
        for said in (desk.get("/api/remote/machines", {}),
                     desk.post("/api/remote/plan", {"config": "study.yml",
                                                    "machine": "box"})):
            assert said["ok"] is False and "hosted GUI" in said["error"]
        assert machine.commands == []

    def test_other_routes_are_not_these(self, machine):
        desk = RemoteDesk(machine.study.parent)
        assert desk.get("/api/status", {}) is None
        assert desk.post("/api/run", {}) is None
        assert desk.get("/api/remote/nothing", {}) is None


STATE = {"system": "1UBQ", "include_phase": ["setup", "simulation"], "output": "ub-run"}


class TestTheConfigBuilder:
    """Run on a machine from the Config Builder: the study is saved as a
    config beside its results folder, as Run on this machine saves it, and
    that file is what is planned and sent."""

    def test_the_study_is_saved_beside_its_results_and_planned(self, served):
        address, machine = served
        root = machine.study.parent
        status, plan = _ask(address, "/api/remote/plan", {"state": STATE, "machine": "box"})
        assert plan["ok"], plan
        assert plan["config"] == "ub-run.yml" and plan["results"] == "ub-run"
        assert plan["fetched_there"] == ["1UBQ"]
        saved = (root / "ub-run.yml").read_text()
        assert "1UBQ" in saved
        assert not _sent(machine)
        # Planned again as it stands: the same file, not another beside it.
        status, again = _ask(address, "/api/remote/plan", {"state": STATE, "machine": "box"})
        assert again["config"] == "ub-run.yml"
        assert sorted(p.name for p in root.glob("ub-run*.yml")) == ["ub-run.yml"]
        # Changed: a file of its own, the first left as it was.
        status, other = _ask(address, "/api/remote/plan", {
            "state": {**STATE, "include_phase": ["setup"]}, "machine": "box"})
        assert other["config"] == "ub-run-2.yml"
        assert (root / "ub-run.yml").read_text() == saved
        status, sent = _ask(address, "/api/remote/send", {"plan": plan["plan"]})
        assert sent["ok"] and sent["job"]["name"] == "ub-run", sent
        travels._until_finished(machine, "ub-run")

    def test_a_results_folder_outside_the_workspace_is_refused(self, served, tmp_path):
        address, machine = served
        status, said = _ask(address, "/api/remote/plan", {
            "state": {**STATE, "output": str(tmp_path / "elsewhere" / "run")},
            "machine": "box"})
        assert said["ok"] is False and "inside the workspace" in said["error"]
        assert not (tmp_path / "elsewhere").exists()
        assert machine.commands == []

    def test_a_window_reaches_the_folder_it_puts_new_studies_in(self, machine, tmp_path):
        """Given the window's runtime, the routes reach the folder it puts
        new studies in as well as the one it opened, and a bare results
        name is where Run on this machine would write it."""
        root = machine.study.parent
        opened = root / "opened"
        opened.mkdir()
        desk = RemoteDesk(DashboardRuntime(workspace_root=opened, exploration_root=root))
        assert desk.roots == [opened.resolve(), root.resolve()]
        assert desk.inside(str(machine.study)) == machine.study.resolve()
        # A relative name is the opened folder's, as a page gives it.
        assert desk.inside("study.yml") == opened.resolve() / "study.yml"
        saved, refused, _ = desk._saved(STATE)
        assert refused is None and saved == root.resolve() / "ub-run.yml"
        assert desk.inside(str(tmp_path.parent)) is None

    def test_a_machine_never_inspected_leaves_no_config_behind(self, served):
        """First review: a plan refused for its machine still saved a
        config, a new one each time the results folder was left blank."""
        address, machine = served
        root = machine.study.parent
        status, said = _ask(address, "/api/remote/plan", {
            "state": {**STATE, "output": ""}, "machine": "gpu-box"})
        assert said["ok"] is False and "gpu-box" in said["error"]
        assert list(root.glob("*.yml")) == [machine.study]

    def test_a_config_saved_for_a_plan_refused_is_taken_back(self, served):
        address, machine = served
        root = machine.study.parent
        (root / "ub-run").mkdir()
        (root / "ub-run" / "kept.txt").write_text("x")
        status, said = _ask(address, "/api/remote/plan", {"state": STATE, "machine": "box"})
        assert said["ok"] is False and said["code"] == "environment.path.exists"
        assert not (root / "ub-run.yml").exists()
        assert (root / "ub-run" / "kept.txt").read_text() == "x"

    def test_a_link_in_the_config_s_place_is_passed_over_unread(self, served, tmp_path):
        address, machine = served
        root = machine.study.parent
        outside = tmp_path / "outside.yml"
        outside.write_text("output: mine\n")
        (root / "ub-run.yml").symlink_to(outside)
        (root / "ub-run-2.yml").write_bytes(b"\xff\xfe")
        status, plan = _ask(address, "/api/remote/plan", {"state": STATE, "machine": "box"})
        assert plan["ok"], plan
        assert plan["config"] == "ub-run-3.yml"
        assert outside.read_text() == "output: mine\n"
        assert (root / "ub-run-2.yml").read_bytes() == b"\xff\xfe"

    def test_a_results_folder_in_the_home_folder_is_refused_before_saving(
            self, machine, tmp_path, monkeypatch):
        home = machine.study.parent
        monkeypatch.setenv("HOME", str(home))
        desk = RemoteDesk(DashboardRuntime(workspace_root=home, exploration_root=home))
        saved, refused, _ = desk._saved(STATE)
        assert saved is None and "home folder" in refused["error"]
        assert not (home / "ub-run.yml").exists()

    def test_a_window_on_a_folder_in_home_reaches_only_that_folder(
            self, served, machine, tmp_path, monkeypatch):
        """Round 1 of the next round: a GUI opened on a folder directly in
        home puts new studies in home, and its routes reached every folder
        in home with them, a config and the files beside it included."""
        home = tmp_path / "h"
        opened = home / "proj"
        opened.mkdir(parents=True)
        notes = home / "Documents" / "notes"
        notes.mkdir(parents=True)
        (notes / "private.pdb").write_text("ATOM secret\n")
        (notes / "c.yml").write_text(machine.study.read_text())
        monkeypatch.setenv("HOME", str(home))
        desk = RemoteDesk(DashboardRuntime(workspace_root=opened, exploration_root=home))
        assert desk.roots == [opened.resolve()]
        assert desk.inside(str(notes / "c.yml")) is None
        said = desk.plan(str(notes / "c.yml"), "box")
        assert said["ok"] is False
        saved, refused, _ = desk._saved(STATE)
        assert saved is None and "home folder" in refused["error"]
        assert not (home / "ub-run.yml").exists()
        assert not _sent(machine)

    def test_a_window_whose_studies_go_in_home_says_so_before_any_plan(
            self, served, tmp_path, monkeypatch):
        """Round 1 of the next round: Plan the send was offered, and only
        once pressed said to open the GUI on a folder of its own, which it
        already was; it now says which folder to open, before anything."""
        home = tmp_path / "h"
        opened = home / "study"
        opened.mkdir(parents=True)
        monkeypatch.setenv("HOME", str(home))
        desk = RemoteDesk(DashboardRuntime(workspace_root=opened, exploration_root=home))
        why = desk.machines()["not_built_here"]
        assert "your home folder" in why and "A config I have" in why
        assert f"`fastmdx gui --output {opened.resolve() / 'first'}`" in why
        saved, refused, _ = desk._saved(STATE)
        assert refused["error"] == why
        # Opened on the home folder itself: the folder above is not called
        # the home folder, and a folder two deeper is named.
        top = RemoteDesk(DashboardRuntime(workspace_root=home, exploration_root=home.parent))
        said = top.machines()["not_built_here"]
        assert "a folder above your home folder" in said and "deeper" not in said
        assert str(home.resolve() / "md" / "first") in said
        # Round 2: opened on a folder outside home whose parent is above it,
        # the folder to open is one inside it, not one moved into home.
        scratch = tmp_path / "scratch"
        scratch.mkdir()
        outside = RemoteDesk(DashboardRuntime(workspace_root=scratch,
                                              exploration_root=Path("/")))
        assert str(scratch.resolve() / "first") in outside.machines()["not_built_here"]
        assert RemoteDesk(DashboardRuntime(workspace_root=opened / "first",
                                           exploration_root=opened)
                          ).machines()["not_built_here"] is None

    def test_a_relative_config_name_checked_elsewhere_is_refused(
            self, served, tmp_path, monkeypatch):
        """Round 1 of the next round: A config I have checks a relative
        name in the folder the GUI was started in, and the plan read it in
        the folder opened, so another file was planned than was checked."""
        address, machine = served
        started = tmp_path / "typed-here"
        started.mkdir()
        (started / "study.yml").write_text(machine.study.read_text())
        monkeypatch.chdir(started)
        status, said = _ask(address, "/api/remote/plan", {"config": "study.yml",
                                                          "machine": "box"})
        assert said["ok"] is False and said["code"] == "remote.config.read_elsewhere"
        assert str(started.resolve()) in said["error"]
        status, full = _ask(address, "/api/remote/plan", {"config": str(machine.study),
                                                          "machine": "box"})
        assert full["ok"], full
        assert not _sent(machine)

    def test_the_refusal_of_a_relative_config_names_each_folder_as_it_is(
            self, served, tmp_path, monkeypatch):
        """Round 3 of the last round: a name in a folder was said as checked
        from that folder, not the one the GUI was started in; a name going
        up a folder was said unnormalised; and one the plan would read
        outside the workspace was not said to be."""
        address, machine = served
        root = machine.study.parent.resolve()
        started = tmp_path / "typed-here"
        (started / "sub").mkdir(parents=True)
        (started / "sub" / "study.yml").write_text(machine.study.read_text())
        (root / "sub").mkdir()
        (root / "sub" / "study.yml").write_text(machine.study.read_text())
        monkeypatch.chdir(started)
        _, deeper = _ask(address, "/api/remote/plan", {"config": "sub/study.yml",
                                                       "machine": "box"})
        error = deeper["error"]
        assert error.startswith(
            f"Checked from the folder the GUI was started in ({started.resolve()}), "
            f"sub/study.yml is {started.resolve() / 'sub' / 'study.yml'}; a plan reads it "
            f"from the folder the GUI was opened on ({root}), as {root / 'sub' / 'study.yml'}.")
        monkeypatch.chdir(started / "sub")
        _, up = _ask(address, "/api/remote/plan", {"config": "../up/study.yml",
                                                   "machine": "box"})
        # Round 2: of a file outside the routes nothing tells whether it is
        # there; a missing one and one there read alike.
        assert f"is {(started / 'up' / 'study.yml').resolve()};" in up["error"]
        # remote-b5's round 1: the copy way was reached only a step later.
        assert up["error"].endswith(
            "Give the full path of a config in the workspace, or, to send a config kept "
            "elsewhere, copy it, with the files it names, into a new folder inside "
            f"{root} and give its full path there.")
        assert (f"a plan reads a relative name only inside the folder the GUI "
                f"was opened on ({root}), and ../up/study.yml leads out of it.") in up["error"]
        assert up["code"] == "remote.config.read_elsewhere"
        # Round 1: started in the folder opened, the two readings are one,
        # and the refusal says only that a plan does not read it.
        monkeypatch.chdir(root)
        _, same = _ask(address, "/api/remote/plan", {"config": "../up/study.yml",
                                                     "machine": "box"})
        assert same["error"].startswith(
            f"A plan reads a relative name only inside the folder the GUI was "
            f"opened on ({root})")
        assert "Checked from" not in same["error"]
        # A name no file can have is refused as that, not as a folder gone.
        _, odd = _ask(address, "/api/remote/plan", {"config": "a\x00b.yml", "machine": "box"})
        assert odd["ok"] is False and "is gone" not in odd["error"]
        assert "is not a config file in" in odd["error"]
        assert not _sent(machine)

    def test_a_relative_config_is_refused_once_the_folder_started_in_is_gone(
            self, served, tmp_path, monkeypatch):
        """Round 3 of the last round: with the folder the GUI was started in
        removed, the name was let through to the plan, which read the
        folder opened, unchecked."""
        address, machine = served

        def gone():
            raise FileNotFoundError("No such file or directory")

        monkeypatch.setattr(Path, "cwd", staticmethod(gone))
        status, said = _ask(address, "/api/remote/plan", {"config": machine.study.name,
                                                          "machine": "box"})
        assert said["ok"] is False and said["code"] == "remote.config.read_elsewhere"
        assert "is gone" in said["error"] and "full path" in said["error"]
        assert not _sent(machine)

    def test_a_relative_config_checked_outside_the_routes_says_how_to_send_it(
            self, served, tmp_path, monkeypatch):
        """Round 3 of the last round: the refusal said to give the full path,
        which is refused too where the file checked is outside the folders
        the routes reach; and a plan's reading of a file not there was named
        as if it were."""
        _, machine = served
        root = machine.study.parent.resolve()
        opened = root / "opened"
        opened.mkdir()
        desk = RemoteDesk(DashboardRuntime(workspace_root=opened, exploration_root=root))
        far = (tmp_path / "far").resolve()
        far.mkdir()
        (far / "study.yml").write_text(machine.study.read_text())
        (far / "top.pdb").write_text((root / "top.pdb").read_text())
        monkeypatch.chdir(far)
        said = desk.plan("study.yml", "box")
        assert said["code"] == "remote.config.read_elsewhere"
        assert "Give its full path" not in said["error"]
        assert (f"The file checked is outside the workspace ({opened}, and {root} where "
                "new studies go), so its full path would be refused too") in said["error"]
        # Round 1: started in its folder, it was told to start the GUI there;
        # round 2: the two ways ran together, and the copy could overwrite.
        # remote-b4's round 3: "here" and the two "it"s read ambiguously.
        assert ("Start the GUI again in the config's folder without `--output`, so it "
                "opens there, or copy the config, with the files it names, into a new "
                f"folder inside {opened} and give its full path there.") in said["error"]
        # A plan's reading inside the routes says a file is not there.
        assert f"opened on ({opened}), where there is no such file." in said["error"]
        # Round 1: the copy given by its bare name was refused as before; by
        # its full path, as said, it plans. So does a GUI opened on its folder.
        (opened / "new").mkdir()
        for name in ("study.yml", "top.pdb"):
            (opened / "new" / name).write_bytes((far / name).read_bytes())
        assert desk.plan("new/study.yml", "box")["code"] == "remote.config.read_elsewhere"
        copied = desk.plan(str(opened / "new" / "study.yml"), "box")
        assert copied["ok"], copied
        on_it = RemoteDesk(DashboardRuntime(workspace_root=far, exploration_root=far.parent))
        assert on_it.plan("study.yml", "box")["ok"]
        # Checked inside the routes, its full path is the way.
        monkeypatch.chdir(root)
        inside = desk.plan("study.yml", "box")
        assert inside["code"] == "remote.config.read_elsewhere"
        assert inside["error"].endswith(
            "Give its full path, so the file planned is the one checked.")
        assert not _sent(machine)

    def test_no_way_is_offered_for_a_file_kept_with_keys_or_settings(
            self, served, tmp_path, monkeypatch):
        """Round 1: a key, a settings file or any file outside the routes was
        told how to be sent (open the GUI on its folder, or copy it in), and
        the answer told such a file from one not there. Round 2: another
        tool's settings under a dot-folder (a token in `~/.config/doctl`),
        and a folder beside `Library`, were still told how."""
        address, machine = served
        keys = tmp_path / "account" / ".ssh"
        keys.mkdir(parents=True)
        (keys / "id_ed25519").write_text("KEY\n")
        (keys / "c.yml").write_text(machine.study.read_text())
        beside = tmp_path / "beside"
        (beside / ".ssh").mkdir(parents=True)
        (beside / "c.yml").write_text(machine.study.read_text())
        tool = tmp_path / "account" / ".config" / "doctl"
        tool.mkdir(parents=True)
        (tool / "config.yaml").write_text("access-token: x\n")
        windows = tmp_path / "account" / "AppData" / "Roaming" / "GitHub CLI"
        windows.mkdir(parents=True)
        (windows / "hosts.yml").write_text("github.com: {oauth_token: x}\n")
        xdg = tmp_path / "xdg"
        (xdg / "glab-cli").mkdir(parents=True)
        (xdg / "glab-cli" / "config.yml").write_text("token: x\n")
        monkeypatch.setenv("XDG_CONFIG_HOME", str(xdg))
        root = machine.study.parent.resolve()
        # Each said as a name that is no config, alike whether it is there.
        for named in (keys / "id_ed25519", keys / "c.yml", beside / "c.yml",
                      tool / "config.yaml", windows / "hosts.yml",
                      xdg / "glab-cli" / "config.yml", Path("/etc/passwd"),
                      tmp_path / "missing.yml"):
            _, said = _ask(address, "/api/remote/plan", {"config": str(named),
                                                         "machine": "box"})
            assert said["error"] == (
                f"{named} is not a config file in the workspace ({root}). To send a config "
                "kept elsewhere, copy it, with the files it names, into a new folder inside "
                f"{root} and give its full path there."), named
        # Round 2: named by a relative path, a file there and one not there
        # read alike, and say no way.
        monkeypatch.chdir(keys)
        for name in ("c.yml", "id_ed25519", "id_rsa"):
            _, typed = _ask(address, "/api/remote/plan", {"config": name, "machine": "box"})
            assert typed["code"] == "remote.config.read_elsewhere"
            assert typed["error"] == (
                f"Checked from the folder the GUI was started in ({keys.resolve()}), {name} "
                f"is {keys.resolve() / name}; a plan reads it from the folder the GUI was "
                f"opened on ({root}), where there is no such file. Give the full path of "
                "a config in the workspace, or, to send a config kept elsewhere, copy it, "
                f"with the files it names, into a new folder inside {root} and give its full "
                "path there."), name
        assert not _sent(machine)

    def test_a_config_kept_by_the_home_folder_is_told_only_to_be_copied(
            self, machine, tmp_path, monkeypatch):
        """Round 2: a config directly in the home folder was told no way to
        be sent. remote-b4's round 3: the way told it apart from a name not
        there, and from a dot-named file of secrets, in a home holding keys;
        and a link in `.ssh` to a config elsewhere was followed."""
        import pwd

        home = tmp_path / "account"
        (home / ".ssh").mkdir(parents=True)
        (home / "other").mkdir()
        (home / "other" / "s.yml").write_text(machine.study.read_text())
        (home / ".ssh" / "l.yml").symlink_to(home / "other" / "s.yml")
        (home / "direct.yml").write_text(machine.study.read_text())
        (home / ".secrets.yml").write_text("token: x\n")
        real = pwd.getpwuid(0)
        account = type(real)((real.pw_name, real.pw_passwd, real.pw_uid, real.pw_gid,
                              real.pw_gecos, str(home), real.pw_shell))
        monkeypatch.setattr(pwd, "getpwuid", lambda uid: account)
        monkeypatch.setenv("HOME", str(home))
        opened = (home / "md" / "first").resolve()
        opened.mkdir(parents=True)
        desk = RemoteDesk(DashboardRuntime(workspace_root=opened, exploration_root=opened.parent))
        for named in (home / "direct.yml", home / ".secrets.yml", home / "missing.yml",
                      home / ".ssh" / "l.yml"):
            said = desk.plan(str(named), "box")
            assert said["error"] == (
                f"{named} is not a config file in {desk._workspace()}. To send a "
                "config kept elsewhere, copy it, with the files it names, into a new folder "
                f"inside {opened} and give its full path there."), named
        # The config a link in `.ssh` leads to is kept apart (remote-b6's
        # round 2); another beside it is told the ways.
        assert "Open the GUI" not in desk.plan(str(home / "other" / "s.yml"), "box")["error"]
        (home / "other" / "t.yml").write_text(machine.study.read_text())
        assert "Open the GUI on the config's folder" in desk.plan(
            str(home / "other" / "t.yml"), "box")["error"]

    def test_a_cloud_folder_or_a_library_in_a_study_is_no_place_of_keys(
            self, machine, tmp_path, monkeypatch):
        """remote-b4's round 3: `Library` anywhere, so a cloud folder under
        the home folder's `~/Library/CloudStorage` and a study's own
        `library` of ligands, was taken for a place of keys, and the config
        told no way to be sent."""
        import pwd

        home = tmp_path / "account"
        cloud = home / "Library" / "CloudStorage" / "Dropbox" / "study"
        cloud.mkdir(parents=True)
        (cloud / "study.yml").write_text(machine.study.read_text())
        (home / "Library" / "Keychains").mkdir()
        (home / "Library" / "Keychains" / "k.yml").write_text("k: x\n")
        project = home / "proj"
        (project / "library").mkdir(parents=True)
        (project / "study.yml").write_text(machine.study.read_text())
        real = pwd.getpwuid(0)
        account = type(real)((real.pw_name, real.pw_passwd, real.pw_uid, real.pw_gid,
                              real.pw_gecos, str(home), real.pw_shell))
        monkeypatch.setattr(pwd, "getpwuid", lambda uid: account)
        monkeypatch.setenv("HOME", str(home))
        opened = (home / "md" / "first").resolve()
        opened.mkdir(parents=True)
        desk = RemoteDesk(DashboardRuntime(workspace_root=opened, exploration_root=opened.parent))
        for config in (cloud / "study.yml", project / "study.yml"):
            assert "Open the GUI on the config's folder" in desk.plan(str(config), "box")[
                "error"], config
        assert "Open the GUI" not in desk.plan(
            str(home / "Library" / "Keychains" / "k.yml"), "box")["error"]

    def test_where_a_link_outside_the_routes_leads_is_not_said(
            self, served, tmp_path, monkeypatch):
        """remote-b4's round 3: a relative name read through a link outside
        the routes said where the link leads, which a name not there does
        not."""
        address, machine = served
        # Two folders down, so the name read from there and from the folder
        # opened are not one path.
        started = tmp_path / "a" / "started"
        started.mkdir(parents=True)
        secret = tmp_path / "keys" / ".config" / "gh"
        secret.mkdir(parents=True)
        (secret / "hosts.yml").write_text("token: x\n")
        other = tmp_path / "other"
        other.mkdir()
        (other / "lnk.yml").symlink_to(secret / "hosts.yml")
        (other / "dangling.yml").symlink_to(tmp_path / "nowhere" / "x.yml")
        monkeypatch.chdir(started)
        answers = set()
        for name in ("lnk.yml", "dangling.yml", "none.yml"):
            _, said = _ask(address, "/api/remote/plan", {"config": f"../../other/{name}",
                                                         "machine": "box"})
            assert said["code"] == "remote.config.read_elsewhere", name
            assert "hosts.yml" not in said["error"] and "nowhere" not in said["error"], name
            answers.add(said["error"].replace(name, "NAME"))
        # The three read alike but for the name typed.
        assert len(answers) == 1, answers

    def test_a_link_crossed_before_going_up_is_not_followed(
            self, served, tmp_path, monkeypatch):
        """remote-b5's round 1: a relative name through a link in `.ssh`
        and then up (`.ssh/lnk/../../x.yml`) was followed, and the answer
        said where the link leads, which a name with no such link does not."""
        address, machine = served
        home = tmp_path / "home"
        started = home / "md"
        started.mkdir(parents=True)
        (home / ".ssh").mkdir()
        far = tmp_path / "far"
        (far / "a" / "b").mkdir(parents=True)
        (far / "proj").mkdir()
        (far / "proj" / "c.yml").write_text(machine.study.read_text())
        (home / ".ssh" / "lnk").symlink_to(far / "a" / "b")
        monkeypatch.chdir(started)
        answers = set()
        for link in ("lnk", "none"):
            _, said = _ask(address, "/api/remote/plan", {
                "config": f"../.ssh/{link}/../../proj/c.yml", "machine": "box"})
            assert said["code"] == "remote.config.read_elsewhere", link
            assert str(far) not in said["error"] and "GUI on the" not in said["error"], link
            answers.add(said["error"].replace(link, "NAME"))
        assert len(answers) == 1, answers

    def test_a_settings_folder_kept_by_a_link_is_kept_apart(
            self, machine, tmp_path, monkeypatch):
        """remote-b5's round 1: a settings folder kept elsewhere by a link
        (`~/.config/gh -> ~/dotfiles/gh`) was offered a way when named where
        it really is."""
        import pwd

        home = tmp_path / "account"
        dotfiles = home / "dotfiles" / "gh"
        dotfiles.mkdir(parents=True)
        (dotfiles / "hosts.yml").write_text("github.com: {oauth_token: x}\n")
        (home / ".config").mkdir()
        (home / ".config" / "gh").symlink_to(dotfiles)
        real = pwd.getpwuid(0)
        account = type(real)((real.pw_name, real.pw_passwd, real.pw_uid, real.pw_gid,
                              real.pw_gecos, str(home), real.pw_shell))
        monkeypatch.setattr(pwd, "getpwuid", lambda uid: account)
        monkeypatch.setenv("HOME", str(home))
        opened = (tmp_path / "work" / "first").resolve()
        opened.mkdir(parents=True)
        desk = RemoteDesk(DashboardRuntime(workspace_root=opened, exploration_root=opened.parent))
        for named in (dotfiles / "hosts.yml", dotfiles / "missing.yml"):
            said = desk.plan(str(named), "box")
            assert said["error"].startswith(f"{named} is not a config file in"), named

    def test_with_no_home_found_every_plan_says_so_alike(
            self, machine, tmp_path, monkeypatch):
        """remote-b5's round 1: with no home folder found, a config inside
        the workspace met an unclassified error from the machine's record,
        and a config in the real home folder read otherwise than a name not
        there."""
        import pwd

        def nowhere(cls=None):
            raise RuntimeError("Could not determine home directory.")

        def no_entry(uid):
            raise KeyError(uid)

        root = machine.study.parent
        desk = RemoteDesk(DashboardRuntime(workspace_root=root, exploration_root=root))
        monkeypatch.setattr(Path, "home", classmethod(nowhere))
        monkeypatch.setattr(pwd, "getpwuid", no_entry)
        for named in (str(machine.study), str(tmp_path / "missing.yml"), "study.yml"):
            said = desk.plan(named, "box")
            assert said == {"ok": False, "error": (
                "No home folder was found, so nothing is sent from any folder. Set "
                "`HOME` to your home folder and start the GUI again."),
                "code": "remote.input.outside"}, named

    def test_only_a_study_s_config_is_told_how_to_be_sent(
            self, machine, tmp_path, monkeypatch):
        """remote-b5's round 2: a file of keys reached by a link from a place
        of keys (`~/.kube -> ~/dotfiles/kube`, a single file linked from
        `~/.config/gh`), or kept by the system, was told how to be sent, as
        a study's config would be; and the way named the folder a link in
        the name leads to."""
        import pwd

        home = tmp_path / "account"
        kube = home / "dotfiles" / "kube"
        kube.mkdir(parents=True)
        (kube / "prod.yaml").write_text(
            "apiVersion: v1\nclusters: []\nusers: [{name: x, user: {token: t}}]\n")
        (home / ".kube").symlink_to(kube)
        dots = home / "dots"
        dots.mkdir()
        (dots / "gh-hosts.yml").write_text("github.com: {oauth_token: x}\n")
        (home / ".config" / "gh").mkdir(parents=True)
        (home / ".config" / "gh" / "hosts.yml").symlink_to(dots / "gh-hosts.yml")
        system = tmp_path / "etc" / "k3s"
        system.mkdir(parents=True)
        (system / "k3s.yaml").write_text("apiVersion: v1\nclusters: []\n")
        far = tmp_path / "far" / "deep"
        far.mkdir(parents=True)
        (far / "c.yml").write_text(machine.study.read_text())
        (home / "proj").mkdir()
        (home / "proj" / "lnk").symlink_to(far)
        real = pwd.getpwuid(0)
        account = type(real)((real.pw_name, real.pw_passwd, real.pw_uid, real.pw_gid,
                              real.pw_gecos, str(home), real.pw_shell))
        monkeypatch.setattr(pwd, "getpwuid", lambda uid: account)
        monkeypatch.setenv("HOME", str(home))
        opened = (tmp_path / "work" / "first").resolve()
        opened.mkdir(parents=True)
        desk = RemoteDesk(DashboardRuntime(workspace_root=opened, exploration_root=opened.parent))
        for named in (kube / "prod.yaml", dots / "gh-hosts.yml", system / "k3s.yaml",
                      kube / "missing.yaml"):
            said = desk.plan(str(named), "box")
            assert said["error"].startswith(f"{named} is not a config file in"), named
        # A study's config behind a link is told the way by its name as given.
        said = desk.plan(str(home / "proj" / "lnk" / "c.yml"), "box")
        assert f"`cd {home / 'proj' / 'lnk'}`" in said["error"]
        assert str(far) not in said["error"]

    def test_a_name_through_a_place_of_keys_into_the_routes_is_not_followed(
            self, served, tmp_path, monkeypatch):
        """remote-b5's round 2: a name through a link in `.ssh` that leads
        into the folders the routes reach was followed, and planned, which
        told what is in `.ssh`."""
        address, machine = served
        root = machine.study.parent.resolve()
        keys = tmp_path / "home" / ".ssh"
        keys.mkdir(parents=True)
        (keys / "lnk").symlink_to(root)
        for named in (str(keys / "lnk" / "study.yml"), str(keys / "none" / "study.yml")):
            _, said = _ask(address, "/api/remote/plan", {"config": named, "machine": "box"})
            assert said["ok"] is False and "is not a config file in" in said["error"], named
        monkeypatch.chdir(tmp_path / "home")
        answers = set()
        for link in ("lnk", "none"):
            _, said = _ask(address, "/api/remote/plan", {"config": f".ssh/{link}/study.yml",
                                                         "machine": "box"})
            assert str(root / "study.yml") not in said["error"].split("; ")[0], link
            answers.add(said["error"].replace(link, "NAME"))
        assert len(answers) == 1, answers
        assert not _sent(machine)

    def test_a_name_read_alike_from_both_folders_is_said_as_named(
            self, served, tmp_path):
        """remote-b5's round 2: where the GUI was started in the folder it
        opened, the refusal spoke of "the file checked", a check it never
        mentions."""
        address, machine = served
        root = machine.study.parent.resolve()
        far = (tmp_path / "far").resolve()
        far.mkdir()
        (far / "c.yml").write_text(machine.study.read_text())
        _, said = _ask(address, "/api/remote/plan", {"config": "../far/c.yml", "machine": "box"})
        assert said["error"].startswith(
            f"A plan reads a relative name only inside the folder the GUI was opened on "
            f"({root}), and ../far/c.yml leads out of it. The file ../far/c.yml names is "
            f"outside the workspace ({root}), so its full path would be refused too: nothing "
            "outside the workspace is sent from here.")
        # remote-b5's round 3: a blank name read as " is not a config file".
        _, blank = _ask(address, "/api/remote/plan", {"config": "  ", "machine": "box"})
        assert blank["error"] == ("No config file was named. Give the full path of the "
                                  "config to send.")

    def test_a_window_on_the_home_folder_is_told_the_copy_way(
            self, machine, tmp_path, monkeypatch):
        """remote-b5's round 2: in a window opened on the home folder the
        copy into a new folder was left out, though such a copy is sent."""
        import pwd

        home = (tmp_path / "account").resolve()
        home.mkdir()
        real = pwd.getpwuid(0)
        account = type(real)((real.pw_name, real.pw_passwd, real.pw_uid, real.pw_gid,
                              real.pw_gecos, str(home), real.pw_shell))
        monkeypatch.setattr(pwd, "getpwuid", lambda uid: account)
        monkeypatch.setenv("HOME", str(home))
        desk = RemoteDesk(DashboardRuntime(workspace_root=home, exploration_root=home.parent))
        said = desk.plan(str(tmp_path / "far" / "none.yml"), "box")
        assert said["error"].endswith(
            "To send a config kept elsewhere, copy it, with the files it names, into a new "
            f"folder inside {home} and give its full path there.")
        above = RemoteDesk(DashboardRuntime(workspace_root=home.parent,
                                            exploration_root=home.parent))
        assert "To send" not in above.plan(str(tmp_path / "far" / "none.yml"), "box")["error"]

    def test_a_config_or_results_in_a_place_of_settings_never_travel(
            self, served, tmp_path, monkeypatch):
        """remote-b5's round 3: in a window opened on the home folder, a
        settings file there (`~/.config/gh/hosts.yml`) was planned and sent
        as a study's config, token and all, and so was a config reached
        through a link to `.ssh`; a results folder in `.ssh` was planned
        too, for a fetch to fill. The same from the command line."""
        import pwd

        from fastmdxplora.refusals import refusal_of
        from fastmdxplora.remote import send

        _, machine = served
        home = (tmp_path / "account").resolve()
        study = machine.study.read_text()
        (home / ".config" / "gh").mkdir(parents=True)
        (home / ".config" / "gh" / "hosts.yml").write_text(study)
        (home / ".config" / "doctl").mkdir()
        (home / ".config" / "doctl" / "config.yaml").write_text(study)
        (home / ".ssh").mkdir()
        (home / ".ssh" / "c.yml").write_text(study)
        (home / "proj").mkdir()
        (home / "proj" / "keys").symlink_to(home / ".ssh")
        (home / "proj" / "c.yml").write_text(study + f"output: {home / '.ssh' / 'run'}\n")
        for name in ("top.pdb",):
            for folder in (home / ".config" / "gh", home / ".config" / "doctl",
                           home / ".ssh", home / "proj"):
                (folder / name).write_text((machine.study.parent / name).read_text())
        real = pwd.getpwuid(0)
        account = type(real)((real.pw_name, real.pw_passwd, real.pw_uid, real.pw_gid,
                              real.pw_gecos, str(home), real.pw_shell))
        monkeypatch.setattr(pwd, "getpwuid", lambda uid: account)
        monkeypatch.setenv("HOME", str(home))
        desk = RemoteDesk(DashboardRuntime(workspace_root=home, exploration_root=home.parent))
        # A place of keys by name (`.config/gh`) or by a link on the way is
        # not looked into: it reads as no config file. Others say the place.
        for named in (home / ".config" / "gh" / "hosts.yml",
                      home / "proj" / "keys" / "c.yml"):
            said = desk.plan(str(named), "box")
            assert said["ok"] is False and "is not a config file in" in said["error"], named
        # remote-b6's round 2: other programs' settings, in the workspace,
        # were said to be no config: the place is named.
        said = desk.plan(str(home / ".config" / "doctl" / "config.yaml"), "box")
        assert said["error"] == (
            f"The config {home / '.config' / 'doctl' / 'config.yaml'} is in "
            f"{home / '.config'}, where keys or settings are kept, so it is not sent. Put "
            "the study's config and its files in a folder of their own.")
        said = desk.plan(str(home / "proj" / "c.yml"), "box")
        assert said["ok"] is False and said["code"] == "remote.input.outside"
        assert said["error"] == (
            f"The results folder {home / '.ssh' / 'run'} is in {home / '.ssh'}, where keys "
            "or settings are kept, so nothing is fetched into it. Set `output` in the "
            "config (or Results folder in the form) to a folder outside such places.")
        assert not _sent(machine)
        # The command line refuses them as well, config and results.
        for config, output in ((home / ".config" / "gh" / "hosts.yml", None),
                               (home / "proj" / "keys" / "c.yml", None),
                               (machine.study, str(home / ".ssh" / "run"))):
            with pytest.raises(Exception) as caught:
                send.prepare(config, "box", output=output)
            found = refusal_of(caught.value)
            assert found.code == "remote.input.outside", config
            assert "where keys or settings are kept" in found.message, config

    def test_settings_kept_by_a_link_or_by_name_never_travel(
            self, served, tmp_path, monkeypatch):
        """remote-b6's round 1: a settings folder kept elsewhere by a link
        (`~/.config/gh -> ~/dotfiles/gh`) was sent from a window on the home
        folder, token and all, and from anywhere when named where it really
        is; a key under `~/.ssh -> ~/dotfiles/ssh` travelled as an input; a
        settings folder named by `XDG_CONFIG_HOME` or kept in `~/Library`,
        and a prepared study in `~/.local/share`, were sent too."""
        import pwd

        from fastmdxplora.refusals import refusal_of
        from fastmdxplora.remote import send
        from fastmdxplora.remote.inputs import gather_inputs, kept_place_of

        _, machine = served
        home = (tmp_path / "account").resolve()
        study = machine.study.read_text()
        pdb = (machine.study.parent / "top.pdb").read_text()
        gh = home / "dotfiles" / "gh"
        gh.mkdir(parents=True)
        (gh / "hosts.yml").write_text(study + "oauth_token: gho_SECRET\n")
        (gh / "top.pdb").write_text(pdb)
        (home / ".config").mkdir()
        (home / ".config" / "gh").symlink_to(gh)
        (home / "dotfiles" / "ssh").mkdir()
        (home / "dotfiles" / "ssh" / "id_ed25519").write_text("KEY\n")
        (home / ".ssh").symlink_to(home / "dotfiles" / "ssh")
        conf = home / "conf" / "doctl"
        conf.mkdir(parents=True)
        (conf / "c.yml").write_text(study)
        support = home / "Library" / "Application Support" / "tool"
        support.mkdir(parents=True)
        (support / "c.yml").write_text(study)
        cloud = home / "Library" / "CloudStorage" / "Dropbox" / "s"
        cloud.mkdir(parents=True)
        real = pwd.getpwuid(0)
        account = type(real)((real.pw_name, real.pw_passwd, real.pw_uid, real.pw_gid,
                              real.pw_gecos, str(home), real.pw_shell))
        monkeypatch.setattr(pwd, "getpwuid", lambda uid: account)
        monkeypatch.setenv("HOME", str(home))
        monkeypatch.setenv("XDG_CONFIG_HOME", str(home / "conf"))
        desk = RemoteDesk(DashboardRuntime(workspace_root=home, exploration_root=home.parent))
        for named in (home / ".config" / "gh" / "hosts.yml", gh / "hosts.yml",
                      conf / "c.yml", support / "c.yml"):
            said = desk.plan(str(named), "box")
            assert said["ok"] is False, named
            with pytest.raises(Exception) as caught:
                send.prepare(named, "box")
            assert refusal_of(caught.value).code == "remote.input.outside", named
        assert not _sent(machine)
        # The person's own files there are not kept apart.
        assert kept_place_of(cloud / "c.yml", a_file=True) is None
        # A key reached through the link is not sent as an input.
        with pytest.raises(Exception) as caught:
            gather_inputs({"systems": [{"system": "top.pdb"}], "keys": "ssh/id_ed25519"},
                          home / "dotfiles")
        assert refusal_of(caught.value).code == "remote.input.outside"
        assert "keys and credentials are kept" in refusal_of(caught.value).message
        # Nor a prepared study kept in a settings folder.
        prepared = home / ".local" / "share" / "tool"
        prepared.mkdir(parents=True)
        (prepared / "manifest.json").write_text('{"phases": ["setup"]}')
        (home / ".local" / "share" / "s.yml").write_text(study)
        with pytest.raises(Exception) as caught:
            gather_inputs({"systems": [{"system": "top.pdb"}],
                           "simulation": {"resume_from": str(prepared)}},
                          home / "proj")
        assert refusal_of(caught.value).code == "remote.input.outside"
        assert "keys and credentials are kept" in refusal_of(caught.value).message

    def test_the_place_named_is_the_folder_and_only_home_s_own_dot_folders(
            self, tmp_path, monkeypatch):
        """remote-b6's round 1: the refusal named the file as the place keys
        are kept, and took the person's own dot folders at any depth
        (`~/projects/.hidden`) for places of settings; a home folder typed
        in other letters on a disk that ignores them was not seen."""
        import pwd

        from fastmdxplora.remote import inputs
        from fastmdxplora.remote.inputs import kept_place_of

        home = (tmp_path / "account").resolve()
        (home / ".ssh" / "s1").mkdir(parents=True)
        real = pwd.getpwuid(0)
        account = type(real)((real.pw_name, real.pw_passwd, real.pw_uid, real.pw_gid,
                              real.pw_gecos, str(home), real.pw_shell))
        monkeypatch.setattr(pwd, "getpwuid", lambda uid: account)
        monkeypatch.setenv("HOME", str(home))
        assert kept_place_of(home / ".ssh" / "s1" / "study.yml", a_file=True) == home / ".ssh"
        assert kept_place_of(home / ".ssh" / "run") == home / ".ssh"
        assert kept_place_of(home / ".config" / "x" / "c.yml", a_file=True) == home / ".config"
        for own in (home / "projects" / ".hidden" / "study.yml",
                    home / "md" / "s2" / ".runs" / "r1" / "c.yml"):
            assert kept_place_of(own, a_file=True) is None, own
        monkeypatch.setattr(inputs.sys, "platform", "darwin")
        typed = Path(str(home).upper()) / ".config" / "x" / "c.yml"
        assert kept_place_of(typed, a_file=True) is not None

    def test_links_into_places_of_keys_are_found_however_they_are_kept(
            self, tmp_path, monkeypatch):
        """remote-b6's round 2: links were looked for only in a few places
        and the search stopped at the first folder it could not read, so a
        key a dotfiles manager linked into `~/.ssh` or `~/.aws`, an app's
        folder in `~/Library/Application Support` linked to Dropbox, a link
        three folders down in `~/.config`, and one listed after an unreadable
        folder, all travelled; a link to the home folder's own parent kept
        everything back, and an XDG variable that is no full path raised or
        kept the folder the send was started in."""
        import os
        import pwd

        from fastmdxplora.refusals import refusal_of
        from fastmdxplora.remote import inputs
        from fastmdxplora.remote.inputs import gather_inputs, kept_place_of

        home = (tmp_path / "account").resolve()
        dots = home / "dotfiles"
        (dots / "ssh").mkdir(parents=True)
        (dots / "ssh" / "id_ed25519").write_text("KEY\n")
        (dots / "aws").mkdir()
        (dots / "aws" / "credentials").write_text("[default]\n")
        (dots / "deep").mkdir()
        (dots / "deep" / "c.yml").write_text("token: x\n")
        (dots / "late").mkdir()
        (dots / "late" / "hosts.yml").write_text("token: x\n")
        (home / ".ssh").mkdir()
        (home / ".ssh" / "id_ed25519").symlink_to(dots / "ssh" / "id_ed25519")
        (home / ".aws").mkdir()
        (home / ".aws" / "credentials").symlink_to(dots / "aws" / "credentials")
        (home / "Dropbox" / "App").mkdir(parents=True)
        (home / "Dropbox" / "App" / "c.yml").write_text("token: x\n")
        support = home / "Library" / "Application Support"
        support.mkdir(parents=True)
        (support / "App").symlink_to(home / "Dropbox" / "App")
        (home / ".config" / "a" / "b").mkdir(parents=True)
        (home / ".config" / "a" / "b" / "c").symlink_to(dots / "deep")
        (home / ".config" / "aaa").mkdir()
        (home / ".config" / "zzz").symlink_to(dots / "late")
        (home / ".foo").symlink_to(home.parent)
        (home / "proj").mkdir()
        unreadable = home / ".config" / "aaa"
        listing = os.scandir

        def scandir(path="."):
            if Path(path) == unreadable:
                raise PermissionError(13, "Permission denied", str(path))
            return listing(path)

        monkeypatch.setattr(inputs.os, "scandir", scandir)
        real = pwd.getpwuid(0)
        account = type(real)((real.pw_name, real.pw_passwd, real.pw_uid, real.pw_gid,
                              real.pw_gecos, str(home), real.pw_shell))
        monkeypatch.setattr(pwd, "getpwuid", lambda uid: account)
        monkeypatch.setenv("HOME", str(home))
        monkeypatch.setenv("XDG_DATA_HOME", ".")
        monkeypatch.setenv("XDG_CONFIG_HOME", "~nosuchuser/x")
        inputs._KEPT_CACHE.clear()
        for kept in (dots / "ssh" / "id_ed25519", dots / "aws" / "credentials",
                     dots / "deep" / "c.yml", dots / "late" / "hosts.yml"):
            assert kept_place_of(kept, a_file=True) is not None, kept
        with pytest.raises(Exception) as caught:
            gather_inputs({"systems": [{"system": "top.pdb"}], "key": "ssh/id_ed25519"}, dots)
        assert refusal_of(caught.value).code == "remote.input.outside"
        # Nothing else is kept back by the link to home's parent, nor by an
        # XDG folder given as `.`.
        monkeypatch.chdir(home / "proj")
        assert kept_place_of(home / "proj" / "c.yml", a_file=True) is None

    def test_the_person_s_own_folders_are_not_kept_by_apps_links(
            self, tmp_path, monkeypatch):
        """remote-b6's round 3: every link found in the home folder's
        settings made its target a place of keys, so the links sandboxed
        apps keep to `~/Documents` and `~/Desktop` (`~/Library/Containers`),
        Wine's to `~/Documents`, and a link in the Trash refused ordinary
        studies on every Mac."""
        import pwd

        from fastmdxplora.remote import inputs
        from fastmdxplora.remote.inputs import kept_place_of

        home = (tmp_path / "account").resolve()
        for own in ("Documents", "Desktop", "md"):
            (home / own).mkdir(parents=True)
        (home / "md" / "project").mkdir()
        data = home / "Library" / "Containers" / "app" / "Data"
        data.mkdir(parents=True)
        (data / "Documents").symlink_to(home / "Documents")
        (data / "Desktop").symlink_to(home / "Desktop")
        wine = home / ".wine" / "drive_c" / "users" / "me"
        wine.mkdir(parents=True)
        (wine / "Documents").symlink_to(home / "Documents")
        trash = home / ".local" / "share" / "Trash" / "files"
        trash.mkdir(parents=True)
        (trash / "data").symlink_to(home / "md" / "project")
        support = home / "Library" / "Application Support"
        support.mkdir(parents=True)
        (support / "App").symlink_to(home / "Documents")
        real = pwd.getpwuid(0)
        account = type(real)((real.pw_name, real.pw_passwd, real.pw_uid, real.pw_gid,
                              real.pw_gecos, str(home), real.pw_shell))
        monkeypatch.setattr(pwd, "getpwuid", lambda uid: account)
        monkeypatch.setenv("HOME", str(home))
        inputs._KEPT_CACHE.clear()
        for study in (home / "Documents" / "lysozyme" / "study.yml",
                      home / "Desktop" / "s" / "study.yml",
                      home / "md" / "project" / "study.yml"):
            assert kept_place_of(study, a_file=True) is None, study

    def test_a_home_whose_name_folds_to_other_letters_is_seen(self, tmp_path, monkeypatch):
        """remote-b6's round 2: a home folder named in letters that fold to
        more (`Stra\u00dfe` for `strasse`) was cut at the wrong place, and its
        own `.local` was missed."""
        import pwd

        from fastmdxplora.remote import inputs
        from fastmdxplora.remote.inputs import kept_place_of

        home = (tmp_path / "strasse").resolve()
        home.mkdir()
        real = pwd.getpwuid(0)
        account = type(real)((real.pw_name, real.pw_passwd, real.pw_uid, real.pw_gid,
                              real.pw_gecos, str(home), real.pw_shell))
        monkeypatch.setattr(pwd, "getpwuid", lambda uid: account)
        monkeypatch.setenv("HOME", str(home))
        monkeypatch.setattr(inputs.sys, "platform", "darwin")
        # Typed with the letter that folds to two.
        typed = home.parent / "Stra\u00dfe" / ".local" / "x" / "c.yml"
        assert kept_place_of(typed, a_file=True) is not None

    def test_only_a_study_s_config_travels_and_each_refusal_says_its_way(
            self, served, tmp_path, monkeypatch):
        """remote-b6's round 3: a file of settings could travel as a config
        wherever no rule of places reached it; the command line's results
        remedy offered a way that fails where `--output` was given; a
        settings folder named with a trailing `/` was not kept; a place one
        cannot read, or one too large to look through, broke or passed the
        send; and the refusal did not say a link led there."""
        import os
        import pwd

        from fastmdxplora.refusals import StudyError, refusal_of
        from fastmdxplora.remote import inputs, send
        from fastmdxplora.remote.inputs import kept_place_of

        _, machine = served
        home = (tmp_path / "account").resolve()
        work = home / "work"
        work.mkdir(parents=True)
        (work / "hosts.yml").write_text("github.com: {oauth_token: x}\n")
        conf = home / "conf"
        (conf / "tool").mkdir(parents=True)
        dots = home / "dotfiles" / "gh"
        dots.mkdir(parents=True)
        (home / ".config").mkdir()
        (home / ".config" / "gh").symlink_to(dots)
        real = pwd.getpwuid(0)
        account = type(real)((real.pw_name, real.pw_passwd, real.pw_uid, real.pw_gid,
                              real.pw_gecos, str(home), real.pw_shell))
        monkeypatch.setattr(pwd, "getpwuid", lambda uid: account)
        monkeypatch.setenv("HOME", str(home))
        monkeypatch.setenv("XDG_CONFIG_HOME", str(conf) + "/")
        inputs._KEPT_CACHE.clear()
        with pytest.raises(StudyError) as caught:
            send.prepare(work / "hosts.yml", "box")
        assert refusal_of(caught.value).code == "config.option.missing_companion"
        assert kept_place_of(conf / "tool" / "c.yml", a_file=True) is not None
        with pytest.raises(StudyError) as caught:
            send.prepare(machine.study, "box", output=str(home / ".ssh" / "run"))
        assert caught.value.args[0].endswith(
            "Give `--output` naming a folder outside such places.")
        said = inputs.kept_place_said(dots / "hosts.yml", kept_place_of(dots / "hosts.yml"))
        assert f"(a link, {home / '.config' / 'gh'}, leads there)" in said
        # A place that cannot be read is passed over; one too large refuses.
        listing = os.scandir

        def scandir(path="."):
            if Path(path) == home / ".config":
                raise PermissionError(13, "Permission denied", str(path))
            return listing(path)

        monkeypatch.setattr(inputs.os, "scandir", scandir)
        inputs._KEPT_CACHE.clear()
        assert kept_place_of(work / "c.yml", a_file=True) is None
        monkeypatch.setattr(inputs.os, "scandir", listing)
        monkeypatch.setattr(inputs, "_LINK_ENTRIES", 3)
        for n in range(5):
            (home / ".config" / f"f{n}").write_text("")
        inputs._KEPT_CACHE.clear()
        with pytest.raises(StudyError) as caught:
            kept_place_of(work / "c.yml", a_file=True)
        assert "nothing is sent" in caught.value.args[0]

    def test_a_key_behind_a_link_from_a_place_of_keys_never_travels(
            self, tmp_path, monkeypatch):
        """remote-b7's round 1: a link from `~/.ssh`, a place of keys in
        `.config` or a home's own dot entry into a folder directly in the
        home folder (`~/.ssh -> ~/keys`) was taken for an app's link into
        `~/Documents`, so a study in `~/keys` sent the key; a continuation
        named by `prepared_from` was told it named no study; and a link
        four folders down in `~/.config/gcloud` was not seen."""
        import pwd

        from fastmdxplora.refusals import StudyError, refusal_of
        from fastmdxplora.remote import inputs
        from fastmdxplora.remote.inputs import gather_inputs, kept_place_of, names_a_study

        home = (tmp_path / "account").resolve()
        for own in ("keys", "gh", "pass", "Documents", "creds"):
            (home / own).mkdir(parents=True)
        (home / "keys" / "id_ed25519").write_text("KEY\n")
        (home / ".ssh").symlink_to(home / "keys")
        (home / ".config" / "gcloud" / "legacy_credentials" / "me").mkdir(parents=True)
        (home / ".config" / "gh").symlink_to(home / "gh")
        (home / ".config" / "gcloud" / "legacy_credentials" / "me" / "adc").symlink_to(
            home / "creds")
        (home / ".password-store").symlink_to(home / "pass")
        # A program's link in `.config` into the person's own folder.
        (home / ".config" / "app").mkdir()
        (home / ".config" / "app" / "docs").symlink_to(home / "Documents")
        real = pwd.getpwuid(0)
        account = type(real)((real.pw_name, real.pw_passwd, real.pw_uid, real.pw_gid,
                              real.pw_gecos, str(home), real.pw_shell))
        monkeypatch.setattr(pwd, "getpwuid", lambda uid: account)
        monkeypatch.setenv("HOME", str(home))
        inputs._KEPT_CACHE.clear()
        for kept in ("keys", "gh", "pass", "creds"):
            assert kept_place_of(home / kept / "study.yml", a_file=True) == home / kept, kept
        assert kept_place_of(home / "Documents" / "s" / "study.yml", a_file=True) is None
        with pytest.raises(StudyError) as caught:
            gather_inputs({"systems": [{"system": "x.pdb"}],
                           "ligand": {"file": "id_ed25519"}}, home / "keys")
        assert refusal_of(caught.value).code == "remote.input.outside"
        said = inputs.kept_place_said(home / "keys" / "study.yml", home / "keys")
        assert f"(a link, {home / '.ssh'}, leads there)" in said
        assert said.endswith(f"Move the study's config and its files out of {home / 'keys'}.")
        assert names_a_study({"simulation": {"prepared_from": "../first"}})

    def test_no_send_is_kept_back_by_programs_own_large_folders(
            self, tmp_path, monkeypatch):
        """remote-b7's round 1: `LOCALAPPDATA` (its `Temp`) and package
        stores such as `~/.local/share/mamba/pkgs` were looked through for
        links, so one too large to look through refused every send; an XDG
        folder set to the home folder was looked through as one; and a
        place too large was looked through again at every check."""
        import os
        import pwd

        from fastmdxplora.refusals import StudyError
        from fastmdxplora.remote import inputs
        from fastmdxplora.remote.inputs import kept_place_of

        home = (tmp_path / "account").resolve()
        local = home / "AppData" / "Local"
        (local / "Temp").mkdir(parents=True)
        pkgs = home / ".local" / "share" / "mamba" / "pkgs"
        pkgs.mkdir(parents=True)
        (home / "Documents" / "st").mkdir(parents=True)
        (home / "Documents" / "ref").symlink_to(home / "Documents" / "st")
        for n in range(5):
            (local / "Temp" / f"t{n}").write_text("")
            (local / f"l{n}").write_text("")
            (pkgs / f"p{n}").mkdir()
            (home / f"h{n}").write_text("")
        real = pwd.getpwuid(0)
        account = type(real)((real.pw_name, real.pw_passwd, real.pw_uid, real.pw_gid,
                              real.pw_gecos, str(home), real.pw_shell))
        monkeypatch.setattr(pwd, "getpwuid", lambda uid: account)
        monkeypatch.setenv("HOME", str(home))
        monkeypatch.setenv("LOCALAPPDATA", str(local))
        monkeypatch.setenv("XDG_CONFIG_HOME", str(home))
        monkeypatch.setattr(inputs, "_LINK_ENTRIES", 3)
        inputs._KEPT_CACHE.clear()
        study = home / "Documents" / "st" / "study.yml"
        assert kept_place_of(study, a_file=True) is None
        # LOCALAPPDATA is still a place of settings by its name.
        assert kept_place_of(local / "x" / "c.yml", a_file=True) is not None
        # A place too large refuses once, and is not looked through again.
        for n in range(5):
            (home / ".config" / f"f{n}").mkdir(parents=True)
        inputs._KEPT_CACHE.clear()
        listing = os.scandir
        looked = []

        def scandir(path="."):
            looked.append(Path(path))
            return listing(path)

        monkeypatch.setattr(inputs.os, "scandir", scandir)
        for _ in range(2):
            with pytest.raises(StudyError) as caught:
                kept_place_of(study, a_file=True)
        assert looked.count(home / ".config") == 1
        assert "holds more than 3 entries" in caught.value.args[0]
        assert "Remove old caches or data you no longer need" in caught.value.args[0]

    def test_a_place_too_large_to_look_through_is_said_in_the_gui(
            self, tmp_path, monkeypatch):
        """remote-b7's round 1: in the GUI the refusal of a place of settings
        too large to look through was read as any error, so a study's config
        outside the window was told to be copied in, which was refused in
        turn for the real reason."""
        import pwd

        from fastmdxplora.remote import inputs

        home = (tmp_path / "account").resolve()
        root = home / "md" / "first"
        root.mkdir(parents=True)
        study = home / "Documents" / "proj" / "study.yml"
        study.parent.mkdir(parents=True)
        study.write_text("systems: [{system: 1L2Y}]\n")
        share = home / ".local" / "share"
        share.mkdir(parents=True)
        for n in range(5):
            (share / f"s{n}").write_text("")
        real = pwd.getpwuid(0)
        account = type(real)((real.pw_name, real.pw_passwd, real.pw_uid, real.pw_gid,
                              real.pw_gecos, str(home), real.pw_shell))
        monkeypatch.setattr(pwd, "getpwuid", lambda uid: account)
        monkeypatch.setenv("HOME", str(home))
        monkeypatch.setattr(inputs, "_LINK_ENTRIES", 3)
        inputs._KEPT_CACHE.clear()
        said = RemoteDesk(root).plan(str(study), "box")
        assert said["ok"] is False
        assert f"{share} holds more than 3 entries" in said["error"], said
        assert "copy the config" not in said["error"]

    def test_a_config_that_cannot_be_read_is_said_without_its_lines(self, served):
        """remote-b6's round 1: a config that could not be read was said with
        the YAML reader's message, which quotes a line of the file."""
        address, machine = served
        bad = machine.study.parent / "bad.yml"
        bad.write_text("token: SECRET_LINE: [\n")
        _, said = _ask(address, "/api/remote/plan", {"config": str(bad), "machine": "box"})
        assert said["ok"] is False and "could not be read" in said["error"]
        assert "SECRET_LINE" not in said["error"]

    def test_a_name_through_a_link_into_a_place_of_keys_is_not_followed(
            self, served, tmp_path, monkeypatch):
        """remote-b5's round 3: a name through a link to `.ssh` with no dot
        in its own parts (`proj/keys/lnk/x.yml`, `keys -> ~/.ssh`) was
        followed into the routes, and told what is in `.ssh`."""
        address, machine = served
        root = machine.study.parent.resolve()
        home = tmp_path / "home"
        (home / ".ssh").mkdir(parents=True)
        (home / ".ssh" / "lnk").symlink_to(root)
        (home / "proj").mkdir()
        (home / "proj" / "keys").symlink_to(home / ".ssh")
        monkeypatch.chdir(home)
        answers = set()
        for link in ("lnk", "none"):
            _, said = _ask(address, "/api/remote/plan", {
                "config": f"proj/keys/{link}/study.yml", "machine": "box"})
            assert said["ok"] is False, link
            answers.add(said["error"].replace(link, "NAME"))
            _, full = _ask(address, "/api/remote/plan", {
                "config": str(home / "proj" / "keys" / link / "study.yml"), "machine": "box"})
            assert "is not a config file in" in full["error"], link
        assert len(answers) == 1, answers
        assert not _sent(machine)

    def test_a_config_is_judged_without_waiting_or_reading_much(self, tmp_path):
        """remote-b5's round 3: a config outside the routes was read whole
        once its size was seen, so a file grown since was read in full and a
        pipe put in its place held the request."""
        import os

        from fastmdxplora.gui.remote_routes import _CONFIG_READ_BYTES, _a_study_config

        pipe = tmp_path / "c.yml"
        os.mkfifo(pipe)
        said: list[bool] = []
        # In a thread, so a read held by the pipe fails the test, not the run.
        reading = threading.Thread(target=lambda: said.append(_a_study_config(pipe)),
                                   daemon=True)
        reading.start()
        reading.join(timeout=5)
        assert said == [False], "held by a pipe"
        big = tmp_path / "big.yml"
        big.write_text("systems: []\n" + "# x\n" * _CONFIG_READ_BYTES)
        assert _a_study_config(big) is False
        small = tmp_path / "small.yml"
        small.write_text("systems: [{system: 1UBQ}]\n")
        assert _a_study_config(small) is True

    def test_a_send_from_the_command_line_with_no_home_found_says_so(
            self, machine, tmp_path, monkeypatch):
        """remote-b4's round 3: with no home folder found, a send was
        refused by the look for your defaults, as an unclassified error,
        before its own refusal could say why."""
        import pwd

        from fastmdxplora.refusals import refusal_of
        from fastmdxplora.remote import send

        def nowhere(cls=None):
            raise RuntimeError("Could not determine home directory.")

        def no_entry(uid):
            raise KeyError(uid)

        monkeypatch.setattr(Path, "home", classmethod(nowhere))
        monkeypatch.setattr(pwd, "getpwuid", no_entry)
        with pytest.raises(Exception) as caught:
            send.prepare(machine.study, "box")
        found = refusal_of(caught.value)
        assert found.code == "remote.input.outside"
        assert "No home folder was found" in found.message

    def test_a_hard_link_to_the_config_planned_is_not_the_one_checked(
            self, served, tmp_path, monkeypatch):
        """Round 1: the same file reached from another folder was taken for
        the one planned, though the files it names by relative paths are
        read from that other folder."""
        address, machine = served
        far = tmp_path / "far"
        far.mkdir()
        (far / "study.yml").hardlink_to(machine.study)
        (far / "top.pdb").write_text("ATOM      other\n")
        monkeypatch.chdir(far)
        _, said = _ask(address, "/api/remote/plan", {"config": "study.yml", "machine": "box"})
        assert said["ok"] is False and said["code"] == "remote.config.read_elsewhere"
        assert not _sent(machine)

    def test_a_command_in_a_note_is_quoted_or_not_shown(self, machine, tmp_path):
        """Round 1: the note's commands named their folders unquoted, so a
        name with a space or `$(...)` broke the command or ran another; and
        a name no quoting holds safely is said in words."""
        from fastmdxplora.gui.remote_routes import _as_typed

        root = machine.study.parent.resolve()
        odd = root / "my $(touch x) runs"
        odd.mkdir()
        desk = RemoteDesk(DashboardRuntime(workspace_root=odd, exploration_root=Path("/")))
        said = desk.machines()["not_built_here"]
        assert f"`fastmdx gui --output '{odd / 'first'}'`" in said
        curly = root / "it\u2019s"
        assert _as_typed(curly) is None
        # Round 2: fish reads a backslash in single quotes as an escape.
        assert _as_typed(root / "a\\';echo x;'") is None
        assert _as_typed(root / "it's") is None
        # remote-b4's round 3: csh reads `!` inside single quotes.
        assert _as_typed(root / "a!b") is None
        curly.mkdir()
        words = RemoteDesk(DashboardRuntime(workspace_root=curly, exploration_root=Path("/")))
        said = words.machines()["not_built_here"]
        assert f"with `fastmdx gui --output` naming a new folder in {curly}." in said

    def test_a_name_the_folder_started_in_cannot_follow_is_refused(
            self, served, tmp_path, monkeypatch):
        """Round 3 of the last round: a name that could not be followed from
        the folder the GUI was started in (a loop of links there) was let
        through to the plan, which read another file."""
        address, machine = served
        root = machine.study.parent
        started = tmp_path / "typed-here"
        started.mkdir()
        (started / "loop").symlink_to(started / "loop")
        (root / "loop").mkdir()
        (root / "loop" / "study.yml").write_text(machine.study.read_text())
        monkeypatch.chdir(started)
        status, said = _ask(address, "/api/remote/plan", {"config": "loop/study.yml",
                                                          "machine": "box"})
        # Refused, whether the loop is said (Python before 3.13) or read
        # through as another file.
        assert said["ok"] is False and said["code"] == "remote.config.read_elsewhere"
        assert not _sent(machine)

    def test_a_config_the_disk_reads_as_the_one_checked_is_planned(
            self, machine, tmp_path, monkeypatch):
        """Round 3 of the last round: where the disk ignores the case of
        letters, the file checked and the one planned were told apart by
        their names' text, and the plan refused."""
        import os

        root = machine.study.parent.resolve()
        desk = RemoteDesk(DashboardRuntime(workspace_root=root, exploration_root=root))
        # A folder of its own: beside the fixture's `study`, a disk that
        # ignores case (macOS's usual one) reads `Study` as that folder,
        # which already exists (CI on b9c66bbd).
        started = tmp_path / "typed-here" / "Study"
        started.mkdir(parents=True)
        (started / "study.yml").write_text("output: other\n")
        monkeypatch.chdir(started)
        real = os.path.samefile
        # The disk reads the folder and the file in other letters as these.
        alike = ({str(started.resolve() / "study.yml"), str(root / "study.yml")},
                 {str(started.resolve()), str(root)})
        monkeypatch.setattr(os.path, "samefile", lambda a, b: (
            {str(a), str(b)} in alike or real(a, b)))
        assert desk._read_elsewhere_here("study.yml") is None

    def test_a_home_folder_that_cannot_be_looked_in_still_opens_the_panel(
            self, machine, tmp_path, monkeypatch):
        """Round 3 of the last round: a home folder whose look raised (no
        permission to search it) broke the machines route, and the panel did
        not open."""
        import pwd

        top = tmp_path / "homes"
        home = top / "account"
        home.mkdir(parents=True)
        real = pwd.getpwuid(0)
        account = type(real)((real.pw_name, real.pw_passwd, real.pw_uid, real.pw_gid,
                              real.pw_gecos, str(home), real.pw_shell))
        monkeypatch.setattr(pwd, "getpwuid", lambda uid: account)
        monkeypatch.setenv("HOME", str(home))
        looked = Path.is_dir

        def is_dir(self, *args, **kwargs):
            if self == home.resolve():
                raise PermissionError(13, "Permission denied", str(self))
            return looked(self, *args, **kwargs)

        monkeypatch.setattr(Path, "is_dir", is_dir)
        desk = RemoteDesk(DashboardRuntime(workspace_root=top, exploration_root=tmp_path))
        said = desk.machines()["not_built_here"]
        assert said.startswith("Studies built in the form would be saved in")
        assert f"`fastmdx gui --output {home.resolve() / 'md' / 'first'}`" in said

    def test_the_settings_under_the_account_s_home_never_travel_whatever_home_says(
            self, tmp_path, monkeypatch):
        """Round 3 of the last round: FastMDXplora's settings were found
        through `HOME` alone, so with `HOME` naming another folder the ones
        under the account's own home were not kept back."""
        import pwd

        from fastmdxplora.remote.inputs import private_in

        home = tmp_path / "account"
        model = home / ".config" / "fastmdxplora" / "model.json"
        model.parent.mkdir(parents=True)
        model.write_text("{}")
        real = pwd.getpwuid(0)
        account = type(real)((real.pw_name, real.pw_passwd, real.pw_uid, real.pw_gid,
                              real.pw_gecos, str(home), real.pw_shell))
        monkeypatch.setattr(pwd, "getpwuid", lambda uid: account)
        monkeypatch.setenv("HOME", str(tmp_path / "elsewhere"))
        monkeypatch.delenv("XDG_CONFIG_HOME", raising=False)
        assert private_in(model) == model
        assert private_in(home / ".config" / "fastmdxplora") is not None

    def test_with_no_home_found_no_config_travels(self, tmp_path, monkeypatch):
        """Round 3 of the last round: with no home folder to be found, a
        config was sent from any folder, though that folder could be the
        home folder."""
        import pwd

        from fastmdxplora.refusals import refusal_of
        from fastmdxplora.remote.inputs import gather_inputs

        def nowhere(cls=None):
            raise RuntimeError("Could not determine home directory.")

        def no_entry(uid):
            raise KeyError(uid)

        monkeypatch.setattr(Path, "home", classmethod(nowhere))
        monkeypatch.setattr(pwd, "getpwuid", no_entry)
        folder = tmp_path / "study"
        folder.mkdir()
        (folder / "top.pdb").write_text("ATOM\n")
        with pytest.raises(Exception) as caught:
            gather_inputs({"systems": [{"system": "top.pdb"}]}, folder)
        found = refusal_of(caught.value)
        assert found.code == "remote.input.outside"
        assert "No home folder was found" in found.message
        assert "Set `HOME` to your home folder" in found.message
        # Round 2: the GUI said nothing of it for a config outside its reach.
        opened = tmp_path / "opened"
        opened.mkdir()
        (folder / "c.yml").write_text("output: c\n")
        desk = RemoteDesk(DashboardRuntime(workspace_root=opened, exploration_root=opened))
        said = desk.plan(str(folder / "c.yml"), "box")
        assert said["error"].endswith("No home folder was found, so nothing is sent from "
                                      "any folder. Set `HOME` to your home folder and start "
                                      "the GUI again.")

    def test_a_results_name_written_in_the_form_is_a_bare_name(self, machine, tmp_path):
        """Round 1 of the next round: the empty Results box was filled with
        the results folder's full path, where Run on this machine reads a
        bare name as the same folder."""
        root = machine.study.parent
        opened = root / "opened"
        opened.mkdir()
        desk = RemoteDesk(DashboardRuntime(workspace_root=opened, exploration_root=root))
        output = root.resolve() / "fastmdxplora_1UBQ_20261010"
        assert desk._results_name(output) == "fastmdxplora_1UBQ_20261010"
        assert desk._results_name(root.resolve() / "deep" / "x") == str(
            root.resolve() / "deep" / "x")

    def test_a_window_with_no_home_found_still_opens_and_sends_no_form_study(
            self, machine, tmp_path, monkeypatch):
        """Round 2 of the next round: with no home folder to be found (no
        HOME, no account entry) the remote routes raised, and the GUI built
        with them did not open."""
        root = machine.study.parent
        opened = root / "opened"
        opened.mkdir()

        def nowhere(cls=None):
            raise RuntimeError("Could not determine home directory.")

        def no_entry(uid):
            raise KeyError(uid)

        import pwd

        monkeypatch.setattr(Path, "home", classmethod(nowhere))
        monkeypatch.setattr(pwd, "getpwuid", no_entry)
        runtime = DashboardRuntime(workspace_root=opened, exploration_root=root)
        with _serving(runtime) as address:
            status, said = _ask(address, "/api/remote/machines")
        assert status == 200 and said["ok"], said
        # Round 3 of the last round: it named no remedy, and said "here" twice.
        assert said["not_built_here"].startswith("No home folder was found")
        assert "Set `HOME` to your home folder and start the GUI again" in said["not_built_here"]
        assert " here" not in said["not_built_here"]
        # Round 3 of the last round: it said "neither `HOME` nor the account
        # names one", untrue where `HOME` is `/`, and offered A config I
        # have, which a send then refuses with no home found.
        assert "neither" not in said["not_built_here"]
        assert "A config I have" not in said["not_built_here"]
        desk = RemoteDesk(runtime)
        assert desk.roots == [opened.resolve()]
        saved, refused, _ = desk._saved(STATE)
        assert saved is None and refused["error"] == said["not_built_here"]

    def test_a_folder_the_disk_reads_as_home_or_a_disk_s_top_is_left_out(
            self, machine, monkeypatch):
        """Round 2 of the next round: home was told by its path's text, so a
        name the disk reads as the home folder (letters in another case, a
        mount of it) and a disk's top were reached as the new-studies
        folder."""
        import os

        from fastmdxplora.remote import inputs

        root = machine.study.parent
        opened = root / "opened"
        opened.mkdir()
        real = os.path.samefile
        monkeypatch.setattr(inputs.os.path, "samefile", lambda a, b: (
            True if (str(a), str(b)) == (str(root.resolve()), str(inputs.Path.home().resolve()))
            else real(a, b)))
        desk = RemoteDesk(DashboardRuntime(workspace_root=opened, exploration_root=root))
        assert desk.roots == [opened.resolve()]
        assert desk.machines()["not_built_here"]
        top = RemoteDesk(DashboardRuntime(workspace_root=opened, exploration_root=Path("/")))
        assert top.roots == [opened.resolve()]
        assert "the top of the disk" in top.machines()["not_built_here"]

    def test_a_wrong_or_empty_home_does_not_let_the_account_s_home_in(
            self, machine, tmp_path, monkeypatch):
        """Round 2 of the last round: home was read from `HOME` alone, so a
        `HOME` naming another folder, or none, let the account's own home
        folder in as the folder new studies go in."""
        import pwd

        root = machine.study.parent
        home = root / "home"
        opened = home / "proj"
        opened.mkdir(parents=True)
        real = pwd.getpwuid(0)
        account = type(real)((real.pw_name, real.pw_passwd, real.pw_uid, real.pw_gid,
                              real.pw_gecos, str(home), real.pw_shell))
        monkeypatch.setattr(pwd, "getpwuid", lambda uid: account)
        for named in (str(tmp_path / "elsewhere"), "", "/", "relative"):
            monkeypatch.setenv("HOME", named)
            desk = RemoteDesk(DashboardRuntime(workspace_root=opened, exploration_root=home))
            assert desk.roots == [opened.resolve()], named
            assert "your home folder" in desk.machines()["not_built_here"]
            saved, refused, _ = desk._saved(STATE)
            assert saved is None and not (home / "ub-run.yml").exists()
            # Round 1: the folder to open is in the account's home, never
            # in a wrong HOME or under the disk's top.
            top = RemoteDesk(DashboardRuntime(workspace_root=home, exploration_root=root))
            said = top.machines()["not_built_here"]
            assert f"`fastmdx gui --output {home.resolve() / 'md' / 'first'}`" in said, named

    def test_the_command_line_keeps_the_account_s_home_out_whatever_home_says(
            self, tmp_path, monkeypatch):
        """Round 2: the account's home is kept out of a send from the command
        line too, where `HOME` names another folder."""
        import pwd

        from fastmdxplora.refusals import refusal_of
        from fastmdxplora.remote.inputs import gather_inputs

        home = tmp_path / "account"
        home.mkdir()
        (home / "top.pdb").write_text("ATOM\n")
        real = pwd.getpwuid(0)
        account = type(real)((real.pw_name, real.pw_passwd, real.pw_uid, real.pw_gid,
                              real.pw_gecos, str(home), real.pw_shell))
        monkeypatch.setattr(pwd, "getpwuid", lambda uid: account)
        monkeypatch.setenv("HOME", str(tmp_path / "elsewhere"))
        with pytest.raises(Exception) as caught:
            gather_inputs({"systems": [{"system": "top.pdb"}]}, home)
        assert refusal_of(caught.value).code == "remote.input.outside"

    def test_a_home_entry_that_cannot_be_read_opens_the_gui(self, machine, tmp_path,
                                                            monkeypatch):
        """Round 2: an account's home that is a loop of links raised out of
        the routes, and the GUI did not open."""
        import pwd

        loop = tmp_path / "loop"
        loop.symlink_to(loop)
        real = pwd.getpwuid(0)
        account = type(real)((real.pw_name, real.pw_passwd, real.pw_uid, real.pw_gid,
                              real.pw_gecos, str(loop), real.pw_shell))
        monkeypatch.setattr(pwd, "getpwuid", lambda uid: account)
        root = machine.study.parent
        desk = RemoteDesk(DashboardRuntime(workspace_root=root, exploration_root=root))
        assert desk.roots == [root.resolve()]

    def test_a_config_is_saved_only_where_the_routes_reach(self, machine, tmp_path):
        """Round 3 of the last round: the form's config was written in the
        window's folder for new studies without asking again whether the
        routes reach it."""
        root = machine.study.parent
        opened = root / "opened"
        opened.mkdir()
        runtime = DashboardRuntime(workspace_root=opened, exploration_root=root)
        desk = RemoteDesk(runtime)
        moved = tmp_path / "moved"
        moved.mkdir()
        runtime.exploration_root = moved
        saved, refused, _ = desk._saved(STATE)
        assert saved is None and "does not send from" in refused["error"]
        assert not list(moved.iterdir())

    def test_a_refusal_names_both_folders_the_routes_reach(self, machine, tmp_path,
                                                           monkeypatch):
        """Round 3 of the last round: refusals called the folder opened the
        workspace, though the routes reach the folder above it as well."""
        root = machine.study.parent.resolve()
        opened = root / "opened"
        opened.mkdir()
        desk = RemoteDesk(DashboardRuntime(workspace_root=opened, exploration_root=root))
        outside = tmp_path / "far" / "c.yml"
        outside.parent.mkdir()
        outside.write_text("output: x\n")
        said = desk.plan(str(outside), "box")
        assert f"the workspace ({opened}, and {root} where new studies go)" in said["error"]
        # And a file that would travel from outside both (a prepared study
        # beside the config's folder, say).
        from types import SimpleNamespace

        from fastmdxplora.remote import api

        far = tmp_path / "far" / "setup"
        planned = SimpleNamespace(inputs=SimpleNamespace(files={"setup": far}))
        monkeypatch.setattr(api, "plan_send", lambda *a, **k: planned)
        sending, refused = desk._planned(opened / "c.yml", "box", opened / "c")
        assert sending is None and refused["code"] == "remote.input.outside"
        assert f"the workspace ({opened}, and {root} where new studies go)" in refused["error"]

    def test_a_name_no_folder_can_have_is_refused(self, machine):
        root = machine.study.parent
        desk = RemoteDesk(DashboardRuntime(workspace_root=root, exploration_root=root))
        saved, refused, _ = desk._saved({**STATE, "output": str(root / "a\0b")})
        assert saved is None and refused["ok"] is False

    def test_the_config_is_saved_where_its_files_are_read_from(self, served):
        """Second review: a results folder named in a subfolder put the
        config there, and a file the study named by a relative path was
        no longer beside it. The window's results name is Run on this
        machine's: a folder beside the others it made."""
        address, machine = served
        root = machine.study.parent
        status, plan = _ask(address, "/api/remote/plan", {
            "state": {**STATE, "output": "newdir/deep"}, "machine": "box"})
        assert plan["ok"], plan
        assert plan["config"] == "newdir_deep.yml" and plan["results"] == "newdir_deep"
        assert (root / "newdir_deep.yml").is_file() and not (root / "newdir").exists()

    def test_a_file_named_relatively_is_sent_only_as_a_run_here_reads_it(
            self, machine, monkeypatch):
        """Third review: a GUI opened on a folder starts Run on this machine
        in the folder above it, so a study naming mine.pdb would have read
        that folder's mine.pdb here and sent the opened folder's. The
        window's routes read it where Run on this machine does; a desk
        given a folder alone refuses the name."""
        from fastmdxplora.gui.server import start_dashboard_session
        from fastmdxplora.remote import api
        from fastmdxplora.remote import send as sending

        parent = machine.study.parent
        opened = parent / "opened"
        opened.mkdir()
        (parent / "mine.pdb").write_text("ATOM  parent\n")
        (opened / "mine.pdb").write_text("ATOM\n")
        monkeypatch.setattr(api, "_link", lambda name, transport: machine.transport())
        monkeypatch.setattr(api, "this_code", lambda: RELEASE)
        monkeypatch.setattr(sending, "this_code", lambda: RELEASE)
        state = {"system": "mine.pdb", "include_phase": ["setup"], "output": "where"}
        session = start_dashboard_session(output=str(opened), host="127.0.0.1", port=0)
        try:
            address = session.url.split("//", 1)[1].rstrip("/")
            status, plan = _ask(address, "/api/remote/plan", {"state": state,
                                                              "machine": "box"})
            assert plan["ok"], plan
            assert [t["from"] for t in plan["travels"]] == [str(parent / "mine.pdb")]
            assert (parent / "where.yml").is_file() and not list(opened.glob("where*.yml"))
            _ask(address, "/api/remote/forget", {"plan": plan["plan"]})
        finally:
            session.server.shutdown()
        saved, said, _ = RemoteDesk(opened)._saved(state)
        assert saved is None and said["code"] == "remote.input.outside"
        assert said["error"].startswith("Run on this machine would read 'mine.pdb' from")
        assert not any(c.startswith("mkdir") for c in machine.commands)

    def test_a_plan_let_go_takes_back_the_config_written_for_it(self, served):
        """Third review: Not now, or a plan taken away as the study changed,
        left its config behind, one more for each blank results name."""
        address, machine = served
        root = machine.study.parent
        status, plan = _ask(address, "/api/remote/plan", {"state": STATE, "machine": "box"})
        assert plan["results_path"] == str((root / "ub-run").resolve())
        status, again = _ask(address, "/api/remote/plan", {"state": STATE, "machine": "box"})
        # Kept while another plan is of it, and taken back with the last.
        assert _ask(address, "/api/remote/forget", {"plan": plan["plan"]})[1]["ok"]
        assert (root / "ub-run.yml").is_file()
        _ask(address, "/api/remote/forget", {"plan": again["plan"]})
        assert not (root / "ub-run.yml").exists()
        status, other = _ask(address, "/api/remote/plan", {
            "state": {**STATE, "output": "other"}, "machine": "box"})
        _ask(address, "/api/remote/forget", {"plan": other["plan"]})
        assert not (root / "other.yml").exists()
        status, sent = _ask(address, "/api/remote/send", {"plan": other["plan"]})
        assert sent["ok"] is False and not _sent(machine)

    def test_a_config_changed_by_hand_or_not_written_for_it_is_kept(self, served):
        address, machine = served
        root = machine.study.parent
        status, plan = _ask(address, "/api/remote/plan", {"state": STATE, "machine": "box"})
        (root / "ub-run.yml").write_text((root / "ub-run.yml").read_text() + "# mine\n")
        _ask(address, "/api/remote/forget", {"plan": plan["plan"]})
        assert (root / "ub-run.yml").read_text().endswith("# mine\n")
        planned = _planned(address)
        _ask(address, "/api/remote/forget", {"plan": planned["plan"]})
        assert machine.study.is_file()
        assert _ask(address, "/api/remote/forget", {"plan": "nothing"})[1] == {"ok": True}

    def test_a_config_another_plan_sent_is_kept_when_this_one_is_let_go(self, served):
        """First review of 1796-1797: two tabs planned one study, the second
        sent it, and Not now in the first deleted the config it was sent as."""
        address, machine = served
        root = machine.study.parent
        status, first = _ask(address, "/api/remote/plan", {"state": STATE, "machine": "box"})
        status, second = _ask(address, "/api/remote/plan", {"state": STATE, "machine": "box"})
        status, sent = _ask(address, "/api/remote/send", {"plan": second["plan"]})
        assert sent["ok"], sent
        _ask(address, "/api/remote/forget", {"plan": first["plan"]})
        assert (root / "ub-run.yml").is_file()
        travels._until_finished(machine, "ub-run")

    def test_a_send_refused_takes_back_the_config_written_for_it(self, served):
        address, machine = served
        root = machine.study.parent
        (root / "mine.pdb").write_text("ATOM  1\n")
        status, plan = _ask(address, "/api/remote/plan", {
            "state": {**STATE, "system": str(root / "mine.pdb")}, "machine": "box"})
        assert plan["ok"], plan
        (root / "mine.pdb").write_text("ATOM  2\n")  # the same size, not the same file
        status, sent = _ask(address, "/api/remote/send", {"plan": plan["plan"]})
        assert sent["code"] == "remote.send.unconfirmed"
        assert not (root / "ub-run.yml").exists() and not _sent(machine)

    def test_a_pdb_id_or_a_word_like_a_folder_is_not_a_file(self, machine, monkeypatch):
        """First review of 1796-1797: a GUI opened on a folder named after
        its protein refused the protein's PDB ID as a file named relatively."""
        from fastmdxplora.gui.server import start_dashboard_session
        from fastmdxplora.remote import api
        from fastmdxplora.remote import send as sending

        parent = machine.study.parent
        opened = parent / "1ubq"
        opened.mkdir()
        (parent / "mine.pdb").write_text("ATOM\n")
        monkeypatch.setattr(api, "_link", lambda name, transport: machine.transport())
        monkeypatch.setattr(api, "this_code", lambda: RELEASE)
        monkeypatch.setattr(sending, "this_code", lambda: RELEASE)
        session = start_dashboard_session(output=str(opened), host="127.0.0.1", port=0)
        try:
            address = session.url.split("//", 1)[1].rstrip("/")
            status, plan = _ask(address, "/api/remote/plan", {"state": {
                "system": "1ubq", "include_phase": ["setup"], "output": "where"},
                "machine": "box"})
            assert plan["ok"], plan
            assert plan["fetched_there"] == ["1ubq"]
            # The study's results, and its config, beside the folder opened,
            # where Run on this machine puts a new study; never inside it.
            assert plan["results_path"] == str((parent / "where").resolve())
            assert not list(opened.iterdir())
        finally:
            session.server.shutdown()
        assert not _sent(machine)


class TestTakingBack:
    def test_only_the_file_written_is_taken_back(self, tmp_path):
        import os

        from fastmdxplora.gui.remote_routes import _sha_of, _taken_back

        written = tmp_path / "a.yml"
        written.write_text("x: 1\n")
        sha = _sha_of(written)
        # Put in its place with the same words: another file, kept.
        other = tmp_path / "b.yml"
        other.write_text("x: 1\n")
        os.replace(other, written)
        _taken_back(written, "0" * 64)
        assert written.is_file()
        sha = _sha_of(written)
        _taken_back(written, sha)
        assert not written.exists()

    def test_a_pipe_in_its_place_is_never_waited_on(self, tmp_path):
        import os

        from fastmdxplora.gui.remote_routes import _sha_of, _taken_back

        pipe = tmp_path / "a.yml"
        os.mkfifo(pipe)
        assert _sha_of(pipe) == ""
        _taken_back(pipe, "0" * 64)
        assert pipe.exists()


class TestTabsAtOnce:
    """Second review of 1796-1798: a plan let go or refused in one tab while
    another tab's send of the same study was still asking its machine took
    back the config that send was sending."""

    @pytest.fixture
    def desk(self, machine, monkeypatch):
        from fastmdxplora.remote import api
        from fastmdxplora.remote import send as sending

        monkeypatch.setattr(api, "_link", lambda name, transport: machine.transport())
        monkeypatch.setattr(api, "this_code", lambda: RELEASE)
        monkeypatch.setattr(sending, "this_code", lambda: RELEASE)
        real = sending.run_here
        monkeypatch.setattr(sending, "run_here", lambda command, runner=None, **more:
                            real(command, runner=machine.local, **more))
        root = machine.study.parent
        return RemoteDesk(DashboardRuntime(workspace_root=root, exploration_root=root)), root

    @staticmethod
    def _held(desk, name):
        """``desk.<name>`` made to wait, once, until let through."""
        asked, through = threading.Event(), threading.Event()
        real = getattr(desk, name)

        def waiting(*args, **kwargs):
            asked.set()
            through.wait(30)
            return real(*args, **kwargs)

        setattr(desk, name, waiting)
        return asked, through, lambda: setattr(desk, name, real)

    def test_a_plan_let_go_while_another_is_sent_leaves_its_config(self, desk, machine):
        desk, root = desk
        first = desk.plan(None, "box", state=STATE)
        second = desk.plan(None, "box", state=STATE)
        asked, through, undo = self._held(desk, "_planned")
        said = {}
        sending = threading.Thread(target=lambda: said.update(desk.send(second["plan"])))
        sending.start()
        assert asked.wait(30)
        undo()
        desk.forget(first["plan"])
        through.set()
        sending.join(60)
        assert said["ok"], said
        assert (root / "ub-run.yml").is_file()
        travels._until_finished(machine, "ub-run")

    def test_a_plan_refused_while_another_of_it_is_kept_leaves_its_config(
            self, desk, machine):
        desk, root = desk
        (root / "ub-run").mkdir()  # in use: the first plan is refused
        (root / "ub-run" / "kept.txt").write_text("x")
        asked, through, undo = self._held(desk, "_plan")
        said = {}
        planning = threading.Thread(
            target=lambda: said.update(desk.plan(None, "box", state=STATE)))
        planning.start()
        assert asked.wait(30)
        undo()
        (root / "ub-run" / "kept.txt").rename(root / "kept.txt")
        kept = desk.plan(None, "box", state=STATE)
        assert kept["ok"], kept
        (root / "kept.txt").rename(root / "ub-run" / "kept.txt")
        through.set()
        planning.join(60)
        assert said["ok"] is False
        assert (root / "ub-run.yml").is_file()
        desk.forget(kept["plan"])
        assert not (root / "ub-run.yml").exists()
        assert not _sent(machine)

    def test_a_plan_dropped_by_age_takes_back_its_config(self, desk, monkeypatch):
        from fastmdxplora.gui import remote_routes

        desk, root = desk
        old = desk.plan(None, "box", state=STATE)
        assert old["ok"] and (root / "ub-run.yml").is_file()
        monkeypatch.setattr(remote_routes, "PLAN_KEPT_S", -1)
        said = desk.send("not-a-plan")
        assert said["code"] == "remote.send.unconfirmed"
        assert not (root / "ub-run.yml").exists()

    def test_a_phase_or_an_analysis_is_a_word_in_a_study_s_folder(self, desk):
        """Second review of 1796-1798: a GUI opened on a finished study's
        folder refused every plan, its `setup` folder read for the phase."""
        desk, root = desk
        (root / "setup").mkdir()
        (root / "rmsd").mkdir()
        said = desk.plan(None, "box", state={
            **STATE, "include_phase": ["setup", "simulation", "analysis"],
            "analysis": {"include": ["rmsd"]}})
        assert said["ok"], said
        assert said["travels"] == []

    def test_a_config_file_plan_refused_takes_back_a_config_let_go_meanwhile(
            self, desk, machine):
        """Round 1 of the next round: the form's plan let go while its
        config was planned as a config file, and that plan then refused,
        left the config behind."""
        desk, root = desk
        first = desk.plan(None, "box", state=STATE)
        asked, through, undo = self._held(desk, "_plan")
        said = {}
        planning = threading.Thread(target=lambda: said.update(
            desk.plan(str(root / "ub-run.yml"), "gpu-box")))
        planning.start()
        assert asked.wait(30)
        undo()
        desk.forget(first["plan"])
        assert (root / "ub-run.yml").is_file()
        through.set()
        planning.join(60)
        assert said["ok"] is False
        assert not (root / "ub-run.yml").exists()

    def test_a_config_file_planned_holds_the_file_the_form_wrote(self, desk, machine):
        """Round 3 of 1796-1799: the form's config planned as a config file
        was taken back by the form's Not now while that plan was made."""
        desk, root = desk
        first = desk.plan(None, "box", state=STATE)
        asked, through, undo = self._held(desk, "_planned")
        said = {}
        planning = threading.Thread(target=lambda: said.update(
            desk.plan(str(root / "ub-run.yml"), "box")))
        planning.start()
        assert asked.wait(30)
        undo()
        desk.forget(first["plan"])
        assert (root / "ub-run.yml").is_file()
        through.set()
        planning.join(60)
        assert said["ok"], said
        desk.forget(said["plan"])
        assert not (root / "ub-run.yml").exists()
