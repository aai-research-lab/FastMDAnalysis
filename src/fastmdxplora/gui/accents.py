"""The interface's five accents, and the one a service hosting the GUI chose.

The colours are in ``gui/static/theme.css``; each is set as ``data-accent``
on a page's root element (Cyan, the default, as none). A person chooses
theirs in the settings popup, kept in the browser; a service sets the
default for its people with ``fastmdx gui --accent``, which is put in
:data:`ACCENT_ENV` so the reports of the studies it starts on this machine
take it too.
"""

from __future__ import annotations

import os

#: The accents, the default first.
ACCENTS = ("cyan", "blue", "violet", "magenta", "graphite")

#: Where `fastmdx gui --accent` leaves the service's accent, for the GUI
#: and for every run it starts.
ACCENT_ENV = "FASTMDX_ACCENT"


def default_accent() -> str:
    """The accent a page opens with before a person chooses one: the
    service's, or Cyan. A name not among :data:`ACCENTS` is read as none."""
    named = os.environ.get(ACCENT_ENV, "").strip().lower()
    return named if named in ACCENTS else ACCENTS[0]
