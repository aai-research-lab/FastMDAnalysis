"""A job that ends as its status is read still says what it left going.

The status script asked whether the job was going, and counted what its run
left going only where it had ended, before it read the exit code. A job that
ended between the two was read failed with its exit code and nothing of the
process its run left working there (CI #749, Ubuntu 3.10: `exit code 1`
where `exit code 1; 1 process it started still works there` was due), and a
run counted as it ended was said as one more process left. What an ended job
left is now counted after its exit code is read, its run gone by then.
"""

from __future__ import annotations

import os
import secrets
import signal
import subprocess
import sys
import time
from pathlib import Path

import pytest

from fastmdxplora.remote.jobs import Job
from fastmdxplora.remote.send import _status_script

pytestmark = pytest.mark.skipif(not Path("/proc/self").is_dir() or sys.platform != "linux",
                                reason="the processes are counted from /proc")

#: The line of the status script that reads the exit code.
READ = "[ -f exit_code ] && echo \"fmdx:exit_code="


def _started(*command: str, env: dict | None = None) -> subprocess.Popen:
    # A session of its own, so it is stopped with all it started.
    return subprocess.Popen(list(command), env=env, stdin=subprocess.DEVNULL,
                            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                            start_new_session=True)


def _seen_in_proc(pid: int, part: str, carried: bytes, within_s: float = 10.0) -> None:
    deadline = time.monotonic() + within_s
    while True:
        try:
            if carried in Path(f"/proc/{pid}/{part}").read_bytes():
                return
        except OSError:
            pass
        assert time.monotonic() < deadline, f"process {pid} never carried {carried!r} in its {part}"
        time.sleep(0.02)


@pytest.fixture()
def going(tmp_path):
    """The job's folder, its script still going (its command line names its
    send's id, as `sh job.sh <id>` does) and a process its run left going
    (carrying the id, as the run's processes do)."""
    folder = tmp_path / "job"
    folder.mkdir()
    # An id of its own: no other process on the machine carries it.
    send_id = secrets.token_hex(8)
    script = _started("sh", "-c", "sleep 60; :", "job.sh", send_id, "x")
    left = _started("sleep", "60", env={**os.environ, "FMDX_SEND_ID": send_id})
    # Each read in /proc as the status script reads it before the job is
    # looked at: on a busy runner a process just started was not yet there
    # with its id, and what it left was counted as nothing (CI, 10-21).
    _seen_in_proc(script.pid, "cmdline", f"job.sh\0{send_id}\0".encode())
    _seen_in_proc(left.pid, "environ", f"FMDX_SEND_ID={send_id}\0".encode())
    yield folder, script, left, send_id
    for process in (script, left):
        try:
            os.killpg(process.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        process.wait()


def _job(folder: Path, script: subprocess.Popen, send_id: str) -> Job:
    return Job(name="trial", machine="box", remote_dir=str(folder), scheduler="process",
               handle=str(script.pid), submitted_at="", code={}, local_output="",
               extra={"send_id": send_id})


def _said(script: str, folder: Path) -> dict[str, list[str]]:
    out = subprocess.run(["sh", "-s"], input=script, capture_output=True, text=True,
                         cwd=folder, timeout=60).stdout
    said: dict[str, list[str]] = {}
    for line in out.splitlines():
        if line.startswith("fmdx:") and "=" in line:
            key, value = line[len("fmdx:"):].split("=", 1)
            said.setdefault(key, []).append(value)
    return said


def test_a_job_ending_between_the_two_says_what_it_left(going):
    folder, script, _left, send_id = going
    asked = _status_script(_job(folder, script, send_id))
    assert asked.count(READ) == 1
    # The job ends just after it was found going, before its exit code is
    # read: its script writes the code.
    said = _said(asked.replace(READ, "echo 1 > exit_code\n" + READ, 1), folder)
    assert said.get("alive") == ["1"]
    assert said.get("exit_code") == ["1"]
    assert said.get("left") == ["1"], said


def test_a_job_already_ended_says_it_once(going):
    folder, script, _left, send_id = going
    script.kill()
    script.wait()
    (folder / "exit_code").write_text("1\n")
    said = _said(_status_script(_job(folder, script, send_id)), folder)
    assert said.get("left") == ["1"], said
    assert "alive" not in said


def test_a_job_still_going_says_nothing_left(going):
    folder, script, _left, send_id = going
    said = _said(_status_script(_job(folder, script, send_id)), folder)
    assert said.get("alive") == ["1"]
    assert "left" not in said and "exit_code" not in said


def test_the_run_ending_as_it_is_counted_is_not_counted_as_left(going):
    """The run carries the send's id too. Counted while it still worked,
    then gone with its exit code written before the count was said, it was
    said as one more process left working there."""
    folder, script, _left, send_id = going
    run = _started("sleep", "60", env={**os.environ, "FMDX_SEND_ID": send_id})
    try:
        asked = _status_script(_job(folder, script, send_id))
        # Just after the scan: the run ends and its exit code is written.
        ended = (f"kill -9 {run.pid}; while ! grep -q '^State:[[:space:]]*Z' "
                 f"/proc/{run.pid}/status 2>/dev/null; do :; done\n"
                 "echo 1 > exit_code\n")
        cut = asked.index('[ "$went" = 1 ] && echo fmdx:alive=1')
        said = _said(asked[:cut] + ended + asked[cut:], folder)
    finally:
        run.kill()
        run.wait()
    assert said.get("exit_code") == ["1"]
    assert said.get("left") == ["1"], said
