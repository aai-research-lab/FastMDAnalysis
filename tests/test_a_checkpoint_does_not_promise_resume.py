"""The checkpoint is written. Nothing here reads it, and three places said
otherwise.

`checkpoint.chk` was described as "for restart / crash recovery" in the
schema, as something that "can be loaded to resume a run from where it left
off" in the runner, and as "Checkpoint: resume the run from here" in the file
browser. There is no `--resume` and no `loadCheckpoint` call anywhere in the
package, so a run that failed at fifteen per cent after fifty minutes left its
work in a file nothing could open -- while three separate places, all of them
somewhere a person goes to find out what is possible, said it could.

Adding the flag is not the fix on its own. The positions and velocities come
back; the biasing state does not. PLUMED does not re-read HILLS unless the
script says RESTART, and none here does, so a resumed metadynamics run would
start again from zero bias inside a well it had already filled and report a
free energy surface that is quietly wrong. A steered pull places its restraint
by absolute step number, so resuming mid-pull puts the anchor somewhere the
protein is not.

This file holds the line at saying only what is true, and records why the
missing feature is more than a missing flag.
"""

from __future__ import annotations

from pathlib import Path

import pytest

SOURCE = Path(__file__).resolve().parents[1] / "src" / "fastmdxplora"


class TestNothingClaimsAResumeThatIsNotThere:
    @pytest.mark.parametrize("module", [
        "config/schema.py",
        "simulation/runner.py",
        "gui/server.py",
    ])
    def test_the_three_places_no_longer_offer_it(self, module: str) -> None:
        text = (SOURCE / module).read_text(encoding="utf-8")
        for claim in ("resume the run from here",
                      "can be loaded to resume a run from where it left"):
            assert claim not in text

    def test_no_resume_flag_is_advertised(self) -> None:
        from fastmdxplora.cli.main import _build_parser

        explore = _build_parser()._subparsers._group_actions[0].choices["explore"]
        options = [s for a in explore._actions for s in a.option_strings]
        assert "--resume" not in options, (
            "if this flag exists now, the notes about RESTART and the moving "
            "restraint need revisiting rather than deleting")

    def test_nothing_loads_a_checkpoint(self) -> None:
        """The claim the other tests rest on. If this ever fails, resume has
        been implemented and the wording should say so again."""
        # A call, not a mention: the schema tells somebody the file can be
        # opened with OpenMM's own `loadCheckpoint`, which is true and useful
        # and is not this software resuming anything.
        offenders = [
            str(path.relative_to(SOURCE))
            for path in SOURCE.rglob("*.py")
            if ".loadCheckpoint(" in path.read_text(encoding="utf-8")
        ]
        assert not offenders, offenders


class TestTheReasonIsRecordedWhereItWillBeFound:
    def test_the_runner_says_why_it_is_not_simple(self) -> None:
        """Next to the code that writes the file, so whoever adds resume
        reads it before starting rather than after."""
        text = (SOURCE / "simulation" / "runner.py").read_text(encoding="utf-8")
        assert "RESTART" in text
        assert "absolute" in text and "step number" in text

    def test_it_names_what_resuming_would_owe_the_manifest(self) -> None:
        """A trajectory assembled from two pieces is not the same object as
        one that ran through, and the analyses read that provenance."""
        text = (SOURCE / "simulation" / "runner.py").read_text(encoding="utf-8")
        assert "manifest" in text


class TestTheCheckpointIsStillWritten:
    def test_the_setting_survives(self) -> None:
        """Withdrawing the promise is not withdrawing the file: it is what a
        crash leaves behind, and OpenMM can open it by hand."""
        from fastmdxplora.config.schema import PHASE_SCHEMAS

        names = [f.name for f in PHASE_SCHEMAS["simulation"].fields]
        assert "checkpoint_interval_steps" in names

    def test_the_reporter_is_still_attached(self) -> None:
        import fastmdxplora.simulation.runner as runner

        assert hasattr(runner, "_attach_checkpoint_reporter")
