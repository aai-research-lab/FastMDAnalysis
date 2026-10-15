"""The tools an AI app reaches other machines with.

Over :mod:`fastmdxplora.remote.api`, so under the rules every interface
keeps: machines only from records made at a terminal, never a prompt, only
the files in a study's own folder sent, a study sent to a workstation only
where its GPUs have room for it, a job asked about at most every 30 s. And two of this server's
own, since what is sent leaves this computer and what comes back can be
large:

- **A study is sent only once the person says so here.** Where the AI app
  cannot put the question to them, nothing is sent: its own approval of the
  call is not enough, as it is for a study started on this machine.
- **Nothing large comes back unasked.** A fetch says each size first; where
  the AI app cannot ask, only one under 100 MB is fetched.

Jobs are those whose results come back into the workspace, so an AI app
given one folder never reads or stops a job sent for another.

The looks, `list_machines` and `remote_status`, are the Agent's too
(:mod:`fastmdxplora.agent.tools`); what is here sends, fetches and stops.
"""

from __future__ import annotations

import time
from pathlib import Path
from typing import Any, NoReturn

from fastmdxplora.mcp.tools import (
    Context,
    ToolError,
    _plan_lines,
    _unused,
    _went_ahead,
    plan_id_of,
)
from fastmdxplora.refusals import refusal_of

__all__ = ["LARGE_BYTES", "start_on_machine"]

#: What a fetch brings without asking, where the AI app cannot ask.
LARGE_BYTES = 100 * 1000 * 1000


def _size(size: int) -> str:
    if size >= 1e9:
        return f"{size / 1e9:.1f} GB"
    return f"{size / 1e6:.1f} MB" if size >= 1_000_000 else f"{size / 1e3:.0f} kB"


def _said_there(exc: BaseException, machine: str = "") -> tuple[str, str]:
    """A refusal from the remote code, with its fix, and its code."""
    found = refusal_of(exc)
    said = found.message
    if found.code == "environment.service.machine_unreachable" and machine:
        # The AI app cannot type a password, and is never given one: a
        # person signs in at a terminal, and the connection kept then is
        # used while it lasts.
        said += (f" A machine that asks for a password or a second factor is "
                 "reached from here only while FastMDXplora's own connection to it "
                 "is open: the person runs `fastmdx remote --machine "
                 f"{machine}` in a terminal and signs in there, which keeps it "
                 "open for ten minutes, then asks again within that time.")
    return said, found.code


def _refuse_there(exc: BaseException, machine: str = "") -> NoReturn:
    """Raise a refusal from the remote code as the tool's own."""
    said, code = _said_there(exc, machine)
    raise ToolError(said, code=code) from None


def _refuse_taken_folder(exc: BaseException, machine: str) -> NoReturn:
    """A refusal of a send whose folder there is taken, with a step an AI
    app can take: it gives no --output and never forces a send."""
    said, code = _said_there(exc, machine)
    if code == "environment.path.exists":
        sentences = said.rstrip().split(". ")
        dropped = ""
        if "--output" in sentences[-1] or "--force-overwrite" in sentences[-1]:
            dropped, sentences = sentences[-1], sentences[:-1]
        # What waiting would settle is kept: a queue that did not answer
        # may in a minute, and a process working there ends. So is how the
        # job holding the folder there is cancelled.
        again = (", or start it again in a minute" if "in a minute" in dropped
                 else ", or start it again once nothing works there" if "works" in said
                 else "")
        cancel_there = dropped.split(", or give another --output")[0].strip()
        kept = [cancel_there] if "`scancel" in cancel_there else []
        said = (". ".join(sentences + kept).rstrip(".") + ". Set `output` in the config to "
                "a new folder name, save it and check it again" + again + ".")
    raise ToolError(said, code=code) from None


def _job_here(ctx: Context, given: str):
    """The job named, if its results come back into the workspace, found
    as the Agent's remote_status finds it."""
    from fastmdxplora.agent.tools import job_here
    from fastmdxplora.mcp.tools import _as_tool_error, _box

    return _as_tool_error(lambda: job_here(_box(ctx), given))


# ---------------------------------------------------------------------------
# Acting
# ---------------------------------------------------------------------------
def start_on_machine(ctx: Context, file: Path, config: dict[str, Any], plan_id: str,
                     machine: str) -> str:
    """`start_study` with `machine`: the checked file sent to run there."""
    from fastmdxplora.remote import api

    requested = str(config.get("output") or file.with_suffix(""))
    where = ctx.workspace.inside(requested)
    if where is None:
        raise ToolError(f"The results folder {requested} is outside the workspace.")
    _unused(ctx, where)
    if _nobody_to_ask(ctx, f"start_study:{file}:{plan_id}:{machine}:"):
        # Refused before the machine is asked anything: a call that cannot
        # end in a send does not reach it.
        _send_unconfirmed(ctx, file, machine, where)
    try:
        sending = api.plan_send(file, machine, output=where)
    except Exception as exc:  # noqa: BLE001 - a refusal, said as one
        _refuse_there(exc, machine)
    if plan_id_of(file) != plan_id:
        # Read again for the plan: the file the person is asked about is
        # the one that was checked.
        raise ToolError(f"{ctx.workspace.shown(file)} changed after it was checked. "
                        "check_study it again and use that plan_id.")
    outside = [source for source in sending.inputs.files.values()
               if ctx.workspace.inside(source) is None]
    if outside:
        raise ToolError(f"{outside[0]} would be sent with the study, and it is outside the "
                        f"workspace ({ctx.workspace.root}); nothing outside it is sent from "
                        "here. Copy it into the config's folder and name it there.",
                        code="remote.input.outside")
    from fastmdxplora.remote.send import room_said, sent_digest, travelling

    if sending.no_room:
        raise ToolError(sending.no_room, code="remote.machine.no_room")
    from fastmdxplora.remote.jobs import job_names

    if sending.job_name in job_names():
        raise ToolError(f"A job called {sending.job_name} was sent from this computer "
                        "before, and a job's name is its results folder's. Set `output` "
                        "in the config to a new folder name, save it and check it again.",
                        code="environment.path.exists")

    shown = ctx.workspace.shown(where)
    prints = travelling(sending)
    travels = [f"  {ctx.workspace.shown(source)} ({_size(prints[name][0])})"
               for name, source in sending.inputs.files.items()]
    message = "\n".join([
        f"Send the study in {ctx.workspace.shown(file)} to {machine} and run it there?",
        *_plan_lines(config),
        f"Runs in: {sending.installation.path[:200]} on {machine}"
        f" ({'SLURM' if sending.scheduler == 'slurm' else 'a detached process'})",
        *(["Sent with it:", *travels] if travels else []),
        *([f"Fetched there from RCSB: {', '.join(sending.inputs.fetched)}"]
          if sending.inputs.fetched else []),
        *room_said(sending, running=_running_said(ctx, sending.running)),
        f"Results come back to {shown} when fetched.",
    ])
    # The answer is to this send: what travels, each file's contents, where
    # it runs.
    sent = sent_digest(sending, prints)
    agreed = _went_ahead(ctx, "send", message,
                         f"start_study:{file}:{plan_id}:{machine}:{where}:{sent}")
    if agreed is None:
        _send_unconfirmed(ctx, file, machine, where)
    if agreed is False:
        return "Not sent: the person did not go ahead."
    if ctx.call is not None and ctx.call.cancelled:
        return "Not sent: the call was cancelled."
    if plan_id_of(file) != plan_id:
        raise ToolError(f"{ctx.workspace.shown(file)} changed while the person was "
                        "asked. check_study it again and show the person that plan.")
    if sent_digest(sending) != sent:
        raise ToolError("A file that travels with the study changed while the person was "
                        "asked, so it is not what they agreed to. Nothing was sent; ask "
                        "them again.", code="remote.send.unconfirmed")
    _unused(ctx, where)
    try:
        job = api.send_planned(sending)
    except Exception as exc:  # noqa: BLE001 - a refusal, said as one
        _refuse_taken_folder(exc, machine)
    how = "SLURM job" if job.scheduler == "slurm" else "process"
    return (f"Sent to {machine} as job {job.name} ({how} {job.handle}); it runs there "
            "on its own, whether or not this AI app stays open. remote_status says how "
            f"it is doing; fetch_study brings its results to {shown} once it has "
            f"finished; cancel_study stops it. Its folder there: {job.remote_dir}.")


def _nobody_to_ask(ctx: Context, about: str) -> bool:
    """Whether this call can neither ask the person nor carries their answer
    to a question about ``about`` (this send of this file to this machine)."""
    call = ctx.call
    if call is None:
        return True
    # An answer counts only where it comes back on a state this server gave
    # out for such a send: another would have the machine asked every call.
    return not call.answers_a_question(about, "send") and not call.can_ask()


def _send_unconfirmed(ctx: Context, file: Path, machine: str, where: Path) -> NoReturn:
    raise ToolError(
        "This AI app cannot ask the person, and a study is sent to another machine "
        "only once they agree to it here. They can send it from a terminal: "
        f"`fastmdx remote send -c {ctx.workspace.shown(file)} --machine {machine} "
        f"--output {ctx.workspace.shown(where)}`.", code="remote.send.unconfirmed")


def _running_said(ctx: Context, running: list[str]) -> str:
    """The studies from here a workstation is running, each named only
    where its results come back here: a job sent for another folder is not
    this AI app's to know of."""
    from fastmdxplora.remote import api

    here = {job.name for job in api.jobs(under=ctx.workspace.root)}
    named = [name for name in running if name in here]
    others = len(running) - len(named)
    if others:
        named.append(f"{others} other stud{'y' if others == 1 else 'ies'}")
    return ", ".join(named)


def _fetch_study(ctx: Context, args: dict[str, Any]) -> str:
    from fastmdxplora.remote import api

    job = _job_here(ctx, args["job"])
    with_trajectory = bool(args.get("with_trajectory"))
    try:
        # The answer kept under 30 s says whether it may have ended; the
        # fetch itself asks the machine again before it copies anything.
        job = api.status(job.name)
        if job.state in ("ready", "running"):
            return (f"{job.name} is still {job.state}"
                    + (f" ({job.detail})" if job.detail else "")
                    + ". It can be fetched once remote_status says done or failed.")
        # Asked once the job has ended, so kept sizes are of what it left.
        sizes = api.fetch_sizes(job.name, max_age_s=api.STATUS_KEPT_S)
    except Exception as exc:  # noqa: BLE001 - a refusal, said as one
        _refuse_there(exc, job.machine)
    if not sizes.run_written:
        tail = job.extra.get("log_tail") or []
        return "\n".join([
            f"{job.name} ended ({job.state}"
            + (f", {job.detail}" if job.detail else "")
            + f") before its run wrote anything on {job.machine}, so there is nothing "
            "to fetch. The end of its log:", *[f"  {line}" for line in tail[-8:]]])
    bringing = sizes.bringing(with_trajectory)
    shown = ctx.workspace.shown(job.local_output)
    staying = (f"Trajectories and checkpoints come too ({sizes.trajectory_files} "
               f"files, {_size(sizes.trajectory)})." if with_trajectory else
               f"Trajectories and checkpoints stay on {job.machine} "
               f"({sizes.trajectory_files} files, {_size(sizes.trajectory)}).")
    message = (f"Fetch {job.name} from {job.machine} into {shown}?\n"
               f"Results: {_size(sizes.results)}. {staying}\n"
               f"In all: {_size(bringing)}.")
    agreed = _went_ahead(ctx, "fetch", message,
                         f"fetch_study:{job.name}:{with_trajectory}:{bringing}")
    if agreed is None and bringing >= LARGE_BYTES:
        raise ToolError(
            f"Fetching {job.name} brings {_size(bringing)}, and this AI app cannot ask "
            "the person first, so it is not fetched from here. They can fetch it from a "
            f"terminal: `fastmdx remote fetch {job.name}`"
            + (" --with-trajectory" if with_trajectory else "") + ".",
            code="remote.fetch.unconfirmed")
    if agreed is False:
        return "Not fetched: the person did not go ahead."
    if ctx.call is not None and ctx.call.cancelled:
        return "Not fetched: the call was cancelled."
    try:
        # No file larger than the whole was said to be, and a margin.
        job, warnings = api.fetch(job.name, with_trajectory=with_trajectory,
                                  most_bytes=bringing + bringing // 10 + 1_000_000)
    except Exception as exc:  # noqa: BLE001 - a refusal, said as one
        _refuse_there(exc, job.machine)
    return "\n".join([f"Fetched {job.name} into {shown} ({_size(bringing)}). "
                      f"read_study reads it.", *warnings])


def _cancel_study(ctx: Context, args: dict[str, Any]) -> str:
    from fastmdxplora.remote import api
    from fastmdxplora.remote.jobs import ABANDONED, FINISHED, load_job
    from fastmdxplora.remote.send import (
        LEFT_STOPPING_S,
        _send_id,
        left_stop_done,
        left_stop_said,
    )

    job = _job_here(ctx, args["job"])
    try:
        # An answer under 30 s old says whether it may still be going; the
        # cancel itself asks the machine again before it stops anything.
        job = api.status(job.name)
    except Exception as exc:  # noqa: BLE001 - a refusal, said as one
        _refuse_there(exc, job.machine)
    if job.state in FINISHED and job.extra.get("queue_silent"):
        return (f"The cluster's queue did not answer, so whether {job.name} is still "
                f"going is not known (last read {job.state}). Nothing was stopped; ask "
                "again in a minute.")
    left = job.state in FINISHED and bool(job.extra.get("left_going"))
    if job.state in FINISHED and not left:
        return f"{job.name} has ended already ({job.state})."
    standing = left_stop_said(job) if left else ""
    if standing == "stopping":
        return (f"{job.name} has ended ({job.state}); what its run left going on "
                f"{job.machine} was asked to stop from here within the last "
                f"{LEFT_STOPPING_S // 60} minutes, and stops at its next frame. Nothing "
                "more was sent.")
    question = (
        f"Stop what {job.name} left going on {job.machine} now? It ended ({job.state}); "
        "processes its run started were asked to stop from here more than "
        f"{LEFT_STOPPING_S // 60} minutes ago and still work there. Stopped now, they "
        "end at once, without a checkpoint; its folder stays."
        if standing == "overdue" else
        f"Stop what {job.name} left going on {job.machine}? It ended ({job.state}), and "
        "processes its run started still work in its folder there; its folder stays."
        if left else
        f"Stop {job.name} on {job.machine}? It stops where it is; its folder there "
        "stays.")
    asked_of = (job.handle, _send_id(job), left, job.extra.get("left_stopped_at"))
    agreed = _went_ahead(ctx, "cancel", question,
                         f"cancel_study:{job.name}:{job.handle}:{_send_id(job)}"
                         + (f":left:{job.extra.get('left_stopped_at') or 0}" if left else ""))
    if agreed is False:
        return "Not stopped: the person did not go ahead."
    if ctx.call is not None and ctx.call.cancelled:
        return "Not stopped: the call was cancelled."
    # What is stopped is what the person was asked about: a job sent again,
    # or a stop asked elsewhere, while they were asked, is asked about anew.
    now = load_job(job.name)
    if (now.handle, _send_id(now), now.state in FINISHED and bool(
            now.extra.get("left_going")), now.extra.get("left_stopped_at")) != asked_of:
        return (f"{job.name} changed while the person was asked (sent again, ended, or "
                "stopped from elsewhere); nothing was stopped. Ask again.")
    asked_at = time.time()
    try:
        job = api.cancel(job.name)
    except Exception as exc:  # noqa: BLE001 - a refusal, said as one
        _refuse_there(exc, job.machine)
    done = left_stop_done(job, asked_at)
    what = (f"what {job.name} left going" if job.state != ABANDONED
            else f"what {job.name}'s run still had going")
    said = {
        "none": f"{job.name} had ended ({job.state}); nothing of its run still works on "
                f"{job.machine}, so nothing was stopped.",
        "stopping": f"{what[0].upper()}{what[1:]} on {job.machine} was asked to stop within "
                    f"the last {LEFT_STOPPING_S // 60} minutes, and stops at its next "
                    "frame; nothing more was sent.",
        "now": f"Stopped {what} on {job.machine} at once, without a checkpoint.",
        "sent": f"Asked {what} on {job.machine} to stop; it stops at its run's next "
                "frame." if job.state != ABANDONED else "",
    }.get(done, "")
    if said:
        return said + f" Its folder there stays at {job.remote_dir}."
    if job.state != ABANDONED:
        return f"{job.name} had ended already ({job.state}); nothing was stopped."
    return (f"Stopped {job.name} on {job.machine}. Its folder there stays at "
            f"{job.remote_dir}.")
