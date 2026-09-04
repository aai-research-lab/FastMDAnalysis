"""Shared fixtures for the suite.

There is deliberately nothing here about optional chemistry backends. An
autouse fixture used to monkeypatch `cli.main._missing_chemistry_backends`
to `lambda: []` for every test, described as neutralising "an early
dependency preflight" that the CLI performs. The CLI does not perform one:
the decision not to was taken in `main()`, the two functions implementing
the other policy were left behind with no caller, and one of them
referenced a name defined nowhere in the tree. The fixture was protecting
the suite from code that could not run.

If a real preflight is added, its test isolation belongs beside it, and it
should be opt-in per test rather than autouse over the whole suite -- an
autouse fixture that disables a guard is indistinguishable from a guard
that does not work.
"""

from __future__ import annotations
