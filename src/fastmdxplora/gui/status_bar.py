"""What the status bar says of this computer and the Agent (statusbar.js).

The strip under every page names the computer the GUI runs on and its GPUs,
and the AI model the Agent asks. The rest of it (the run going here, the
jobs on other machines, the open study's problems) the page has already from
the routes it reads for them.

Read from what is here: the GPUs as ``nvidia-smi`` gives them
(:func:`fastmdxplora.gpu_here.room_here`), kept a minute so the strip's polls
do not ask it each time, and the fastest OpenMM platform found, kept for the
life of the server; the AI model as Settings stored it, never its key. On
loopback only (the server's gate): a computer's name and its GPUs are not
said to the network.
"""

from __future__ import annotations

import os
import platform
import socket
import threading
import time
from typing import Any

__all__ = ["status_bar", "forget"]

#: How long the GPUs as read are kept.
GPUS_KEPT_S = 60.0

_lock = threading.Lock()
_gpus: tuple[float, list[dict[str, Any]], bool] | None = None
_platform: str | None = None


def _read_gpus() -> tuple[list[dict[str, Any]], bool]:
    """Each GPU's name and memory, and whether they could be read here."""
    from fastmdxplora.gpu_here import nvidia_smi_here, room_here

    if not nvidia_smi_here():
        return [], False
    room = room_here()
    if room is None:
        return [], True
    return [{"index": gpu.index, "name": gpu.name, "total_mb": gpu.total_mb,
             "free_mb": gpu.free_mb, "busy_pct": gpu.busy_pct} for gpu in room.gpus], True


def _gpus_now() -> tuple[list[dict[str, Any]], bool]:
    global _gpus
    with _lock:
        kept = _gpus
    if kept is not None and time.monotonic() - kept[0] < GPUS_KEPT_S:
        return kept[1], kept[2]
    try:
        gpus, read = _read_gpus()
    except Exception:  # noqa: BLE001 - a strip that cannot say its GPUs says the rest
        gpus, read = [], False
    with _lock:
        _gpus = (time.monotonic(), gpus, read)
    return gpus, read


def forget() -> None:
    """What is kept, forgotten: the GPUs and the platform read again."""
    global _gpus, _platform
    with _lock:
        _gpus = None
    _platform = None


def _fastest_platform() -> str:
    """The OpenMM platform a study here would run on, fastest first; "" where
    OpenMM is not installed."""
    global _platform
    if _platform is not None:
        return _platform
    try:
        import openmm

        names = {openmm.Platform.getPlatform(i).getName()
                 for i in range(openmm.Platform.getNumPlatforms())}
    except Exception:  # noqa: BLE001 - no OpenMM, no platform to name
        names = set()
    _platform = next((name for name in ("CUDA", "HIP", "OpenCL", "CPU") if name in names), "")
    return _platform


def _machine() -> dict[str, Any]:
    gpus, read = _gpus_now()
    return {
        "host": socket.gethostname().split(".")[0],
        "system": platform.system(),
        "cpus": os.cpu_count(),
        "gpus": gpus,
        # Whether nvidia-smi was there to read the GPUs: none read where it
        # is is "no GPU here", where it is not "not read".
        "gpus_read": read,
        "platform": _fastest_platform(),
    }


def _agent() -> dict[str, Any] | None:
    from fastmdxplora.agent.models import PROVIDERS, load_choice

    choice = load_choice()
    if choice is None:
        return None
    return {"provider": choice.provider,
            "provider_label": str(PROVIDERS.get(choice.provider, {}).get("label") or choice.provider),
            "model": choice.model}


def status_bar() -> dict[str, Any]:
    """The computer here and the Agent's AI model, for the status bar."""
    try:
        agent = _agent()
    except Exception:  # noqa: BLE001 - the strip says the rest
        agent = None
    return {"ok": True, "machine": _machine(), "agent": agent}
