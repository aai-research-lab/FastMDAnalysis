"""A study stopped again before its extension's first checkpoint is carried on.

Killed in production, a study is carried on in a segment of its own from
its last checkpoint. Killed again before that segment wrote a checkpoint,
`fastmdx resume` read the segment as the place to go on from, found no
checkpoint there, took the study for one whose production had not begun
and ran it again from its start: the production already run was
overwritten. The segment is now set aside, kept under a name nothing reads
as a segment, and the study is carried on from the piece before it, as long
as it was asked to run. A first piece stopped before its first checkpoint is
still run again from its start; an extended study nothing can carry on is
refused, not run again over its pieces, however often it is asked.
"""

from __future__ import annotations

import json
import shutil
import tempfile
import unittest
from pathlib import Path

import pytest

try:
    import mdtraj  # noqa: F401
    import openmm  # noqa: F401
    import pdbfixer  # noqa: F401
    HAS_BACKENDS = True
except ImportError:  # pragma: no cover - the backends are optional
    HAS_BACKENDS = False

from tests.test_a_carried_on_study_reads_its_energy_whole import _config
from tests.test_a_real_study_runs_end_to_end import TRI_ALANINE


@unittest.skipUnless(HAS_BACKENDS, "OpenMM, PDBFixer and MDTraj are needed")
class TestAStudyStoppedBeforeItsExtensionSCheckpoint(unittest.TestCase):
    """Killed 250 steps into production (its checkpoint at 200), carried on,
    and killed again 40 steps into the segment, before its checkpoint."""

    @classmethod
    def setUpClass(cls):
        from unittest import mock

        pytest.importorskip("openmm.app")
        import openmm.app

        from fastmdxplora import FastMDXplora
        from fastmdxplora.simulation.resume import resume_study

        root = _folder()
        cls.addClassCleanup(shutil.rmtree, root, True)
        pdb = root / "tri-ala.pdb"
        pdb.write_text(TRI_ALANINE)
        real_step = openmm.app.Simulation.step
        taken = {"steps": 0, "limit": 200 + 250}

        def step_until_killed(simulation, steps):
            if taken["steps"] >= taken["limit"]:
                raise RuntimeError("killed")
            taken["steps"] += int(steps)
            return real_step(simulation, steps)

        cls.study = root / "study"
        with mock.patch.object(openmm.app.Simulation, "step", step_until_killed):
            FastMDXplora(config_data=_config(pdb, telemetry_interval=25),
                         output_dir=str(cls.study)).explore()
            # A kill leaves no Manifest.
            (cls.study / "manifest.json").unlink()
            taken.update(steps=0, limit=40)
            cls.stopped = resume_study(cls.study)
        (cls.study / "manifest.json").unlink(missing_ok=True)
        cls.first = cls.study / "simulation" / "production.dcd"
        cls.first_bytes = cls.first.read_bytes()
        cls.segment_had_checkpoint = (
            cls.study / "segment-001" / "simulation" / "checkpoint.chk").is_file()
        cls.answer = resume_study(cls.study)

    def test_the_segment_stopped_before_its_checkpoint(self):
        self.assertFalse(self.stopped["ok"])
        self.assertFalse(self.segment_had_checkpoint)

    def test_it_is_carried_on_not_run_again_from_its_start(self):
        self.assertEqual((self.answer["ok"], self.answer["did"]), (True, "continued"),
                         self.answer.get("error"))
        # The first piece's production is as it was written.
        self.assertEqual(self.first.read_bytes(), self.first_bytes)

    def test_the_segment_is_kept_aside_and_run_again(self):
        aside = self.study / "segment-001-stopped-before-its-checkpoint"
        self.assertEqual(self.answer["set_aside"], str(aside))
        self.assertTrue((aside / "simulation").is_dir())
        self.assertTrue((self.study / "segment-001" / "simulation" / "checkpoint.chk").is_file())
        joined = json.loads((self.study / "joined" / "joined.json").read_text())
        self.assertEqual(joined["segments"], [0, 1])
        self.assertEqual(joined["trimmed"], {"0": 4})

    def test_its_analysis_reads_the_whole_run_once(self):
        record = json.loads((self.study / "analysis" / "thermodynamics" / "options.json")
                            .read_text())
        # 200 steps of the first piece and the 100 the segment ran again, a
        # row every 25: 12.
        self.assertEqual(record["findings"]["thermodynamics"]["samples"], 12)


def _folder() -> Path:
    """A new temporary folder, resolved as the resume resolves a study: on
    macOS the temporary folder is under /var, a link to /private/var."""
    return Path(tempfile.mkdtemp()).resolve()


def _killed(study, pdb, *limits):
    """The study run, killed after each count of steps in turn: the first
    as it is explored, each next one as it is carried on."""
    pytest.importorskip("openmm.app")
    from unittest import mock

    import openmm.app

    from fastmdxplora import FastMDXplora
    from fastmdxplora.simulation.resume import resume_study

    real_step = openmm.app.Simulation.step
    taken = {"steps": 0, "limit": limits[0]}

    def step_until_killed(simulation, steps):
        if taken["steps"] >= taken["limit"]:
            raise RuntimeError("killed")
        taken["steps"] += int(steps)
        return real_step(simulation, steps)

    with mock.patch.object(openmm.app.Simulation, "step", step_until_killed):
        FastMDXplora(config_data=_config(pdb, telemetry_interval=25),
                     output_dir=str(study)).explore()
        for limit in limits[1:]:
            (study / "manifest.json").unlink(missing_ok=True)
            taken.update(steps=0, limit=limit)
            resume_study(study)
    (study / "manifest.json").unlink(missing_ok=True)


@unittest.skipUnless(HAS_BACKENDS, "OpenMM, PDBFixer and MDTraj are needed")
class TestAFirstPieceBeforeItsCheckpointIsRunAgain(unittest.TestCase):
    """Killed 75 steps into production, a frame written and no checkpoint
    yet: run again from its start, as the fix offered for it says, not
    refused."""

    @classmethod
    def setUpClass(cls):
        pytest.importorskip("openmm.app")
        from fastmdxplora.simulation.resume import resume_study

        root = _folder()
        cls.addClassCleanup(shutil.rmtree, root, True)
        pdb = root / "tri-ala.pdb"
        pdb.write_text(TRI_ALANINE)
        cls.study = root / "study"
        _killed(cls.study, pdb, 200 + 75)
        simulation = cls.study / "simulation"
        cls.had = ((simulation / "production.dcd").is_file(),
                   (simulation / "checkpoint.chk").is_file())
        cls.answer = resume_study(cls.study)

    def test_it_stopped_with_a_frame_and_no_checkpoint(self):
        self.assertEqual(self.had, (True, False))

    def test_it_is_run_again_from_its_start(self):
        self.assertEqual((self.answer["ok"], self.answer["did"]), (True, "restarted"),
                         self.answer.get("error"))
        self.assertNotIn("set_aside", self.answer)


@unittest.skipUnless(HAS_BACKENDS, "OpenMM, PDBFixer and MDTraj are needed")
class TestAnExtendedStudyWithNoCheckpointIsNotRunOver(unittest.TestCase):
    """Extended, its extension stopped before its checkpoint and the first
    piece's checkpoint gone: nothing can carry it on, and run again from its
    start it would overwrite both pieces. Refused, the segment still set
    aside and said so; and refused again when asked again, as a service that
    resumes every job it restarts asks, not taken then for a study never
    extended and run over every piece."""

    @classmethod
    def setUpClass(cls):
        pytest.importorskip("openmm.app")
        from fastmdxplora.simulation.resume import resume_study

        root = _folder()
        cls.addClassCleanup(shutil.rmtree, root, True)
        pdb = root / "tri-ala.pdb"
        pdb.write_text(TRI_ALANINE)
        cls.study = root / "study"
        _killed(cls.study, pdb, 200 + 250, 40)
        for name in ("checkpoint.chk", "checkpoint.chk.json", "checkpoint.chk.sha256"):
            (cls.study / "simulation" / name).unlink(missing_ok=True)
        cls.first = cls.study / "simulation" / "production.dcd"
        cls.first_bytes = cls.first.read_bytes()
        cls.answer = resume_study(cls.study)
        cls.again = resume_study(cls.study)

    def test_it_is_refused_and_its_frames_kept(self):
        self.assertEqual((self.answer["ok"], self.answer["did"]), (False, "nothing"))
        self.assertIn("would overwrite the frames of every piece", self.answer["error"])
        self.assertEqual(self.first.read_bytes(), self.first_bytes)

    def test_asked_again_it_is_refused_again(self):
        self.assertEqual((self.again["ok"], self.again["did"]), (False, "nothing"),
                         self.again.get("detail"))
        self.assertEqual(self.first.read_bytes(), self.first_bytes)
        aside = self.study / "segment-001-stopped-before-its-checkpoint"
        self.assertTrue(any((aside / "simulation").glob("*.dcd")))

    def test_the_segment_set_aside_is_said(self):
        aside = self.study / "segment-001-stopped-before-its-checkpoint"
        self.assertEqual(self.answer.get("set_aside"), str(aside))
        self.assertTrue((aside / "simulation").is_dir())


@unittest.skipUnless(HAS_BACKENDS, "OpenMM, PDBFixer and MDTraj are needed")
class TestAnExtensionPastThePlanKeepsItsLength(unittest.TestCase):
    """Killed in production, then extended past its plan (0.0008 ns more
    than its checkpoint, the plan being 0.0006 ns in all) and killed again
    before the extension's checkpoint: carried on to the length asked for,
    not cut back to the plan without a word."""

    @classmethod
    def setUpClass(cls):
        from unittest import mock

        pytest.importorskip("openmm.app")
        import openmm.app

        from fastmdxplora import FastMDXplora
        from fastmdxplora.simulation.resume import (
            extend_study,
            production_done_ns,
            resume_study,
        )

        root = _folder()
        cls.addClassCleanup(shutil.rmtree, root, True)
        pdb = root / "tri-ala.pdb"
        pdb.write_text(TRI_ALANINE)
        real_step = openmm.app.Simulation.step
        taken = {"steps": 0, "limit": 200 + 250}

        def step_until_killed(simulation, steps):
            if taken["steps"] >= taken["limit"]:
                raise RuntimeError("killed")
            taken["steps"] += int(steps)
            return real_step(simulation, steps)

        cls.study = root / "study"
        with mock.patch.object(openmm.app.Simulation, "step", step_until_killed):
            FastMDXplora(config_data=_config(pdb, telemetry_interval=25),
                         output_dir=str(cls.study)).explore()
            (cls.study / "manifest.json").unlink()
            taken.update(steps=0, limit=40)
            extend_study(cls.study, more_ns=0.0008)
        (cls.study / "manifest.json").unlink(missing_ok=True)
        cls.answer = resume_study(cls.study)
        cls.done = production_done_ns(cls.study)

    def test_it_runs_to_the_length_asked_for(self):
        self.assertEqual((self.answer["ok"], self.answer["did"]), (True, "continued"),
                         self.answer.get("error"))
        # 200 steps of 2 fs to the checkpoint, and the 0.0008 ns asked.
        self.assertAlmostEqual(self.done, 0.0012, places=6)


@unittest.skipUnless(HAS_BACKENDS, "OpenMM, PDBFixer and MDTraj are needed")
class TestASecondExtensionPastAFinishedPlanKeepsItsLength(unittest.TestCase):
    """Extended past its plan by 0.0008 ns to 0.0012 ns, the plan being
    0.0006 ns; extended by 0.0008 ns again and killed before that
    segment's checkpoint. Nothing is left of the plan, and the extension
    was dropped with it: the study was said whole, its join refused, and
    so at every resume after. Carried on to the length asked, and found
    still by a resume after one that stopped short of carrying it on, the
    study moved to another folder between them."""

    @classmethod
    def setUpClass(cls):
        from unittest import mock

        pytest.importorskip("openmm.app")
        import openmm.app

        from fastmdxplora import FastMDXplora
        from fastmdxplora.simulation import resume

        root = _folder()
        cls.addClassCleanup(shutil.rmtree, root, True)
        pdb = root / "tri-ala.pdb"
        pdb.write_text(TRI_ALANINE)
        real_step = openmm.app.Simulation.step
        taken = {"steps": 0, "limit": 200 + 250}

        def step_until_killed(simulation, steps):
            if taken["steps"] >= taken["limit"]:
                raise RuntimeError("killed")
            taken["steps"] += int(steps)
            return real_step(simulation, steps)

        cls.study = root / "study"
        with mock.patch.object(openmm.app.Simulation, "step", step_until_killed):
            FastMDXplora(config_data=_config(pdb, telemetry_interval=25),
                         output_dir=str(cls.study)).explore()
            (cls.study / "manifest.json").unlink()
            taken.update(steps=0, limit=10 ** 9)
            resume.extend_study(cls.study, more_ns=0.0008)
            taken.update(steps=0, limit=40)
            resume.extend_study(cls.study, more_ns=0.0008)
        # The Manifest the first extension's join wrote is left, as a kill
        # leaves it: every phase done, so the study read as finished.
        cls.second_had_checkpoint = (
            cls.study / "segment-002" / "simulation" / "checkpoint.chk").is_file()
        # A resume that set the segment aside and stopped there, its own
        # run refused: the next one still finds the length asked.
        with mock.patch.object(resume, "extend_study",
                               lambda *a, **k: {"ok": False, "error": "stopped"}):
            cls.stopped = resume.resume_study(cls.study)
        cls.study = cls.study.rename(root / "moved")
        cls.answer = resume.resume_study(cls.study)
        cls.done = resume.production_done_ns(cls.study)

    def test_the_second_segment_stopped_before_its_checkpoint(self):
        self.assertFalse(self.second_had_checkpoint)
        self.assertEqual(self.stopped["did"], "continued")

    def test_it_runs_to_the_length_asked_for(self):
        self.assertEqual((self.answer["ok"], self.answer["did"]), (True, "continued"),
                         self.answer.get("error"))
        # 200 steps of 2 fs to the first checkpoint, and 0.0008 ns twice.
        self.assertAlmostEqual(self.done, 0.0020, places=6)


def _stepping(limit):
    """OpenMM's steps, stopped as a kill stops them once ``limit`` are
    taken; ``limit`` is a dict, changed between the runs."""
    openmm_app = pytest.importorskip("openmm.app")

    real_step = openmm_app.Simulation.step

    def step_until_killed(simulation, steps):
        if limit["steps"] >= limit["limit"]:
            raise RuntimeError("killed")
        limit["steps"] += int(steps)
        return real_step(simulation, steps)

    return step_until_killed


@unittest.skipUnless(HAS_BACKENDS, "OpenMM, PDBFixer and MDTraj are needed")
class TestAFinishedStudyExtendedAndStoppedIsCarriedOn(unittest.TestCase):
    """A study that finished (0.0006 ns), extended by 0.0008 ns and killed
    before the extension's first checkpoint, its Manifest left as the
    finished study wrote it. Read as finished: "nothing to carry on", the
    segment left in place, and extending it again refused for it. Carried on
    to the length asked, then read as finished once it is."""

    @classmethod
    def setUpClass(cls):
        from unittest import mock

        pytest.importorskip("openmm.app")
        import openmm.app

        from fastmdxplora import FastMDXplora
        from fastmdxplora.simulation.resume import (
            extend_study,
            production_done_ns,
            resume_study,
        )

        root = Path(tempfile.mkdtemp())
        cls.addClassCleanup(shutil.rmtree, root, True)
        pdb = root / "tri-ala.pdb"
        pdb.write_text(TRI_ALANINE)
        cls.study = root / "study"
        taken = {"steps": 0, "limit": 10 ** 9}
        with mock.patch.object(openmm.app.Simulation, "step", _stepping(taken)):
            FastMDXplora(config_data=_config(pdb, telemetry_interval=25),
                         output_dir=str(cls.study)).explore()
            taken.update(steps=0, limit=40)
            extend_study(cls.study, more_ns=0.0008)
        cls.had_manifest = (cls.study / "manifest.json").is_file()
        cls.answer = resume_study(cls.study)
        cls.done = production_done_ns(cls.study)
        cls.again = resume_study(cls.study)

    def test_it_is_carried_on_to_the_length_asked(self):
        self.assertTrue(self.had_manifest)
        self.assertEqual((self.answer["ok"], self.answer["did"]), (True, "continued"),
                         self.answer.get("error") or self.answer.get("detail"))
        aside = self.study / "segment-001-stopped-before-its-checkpoint"
        self.assertEqual(self.answer.get("set_aside"), str(aside))
        # 300 steps of 2 fs, and the 0.0008 ns asked.
        self.assertAlmostEqual(self.done, 0.0014, places=6)

    def test_once_carried_on_it_is_finished(self):
        self.assertEqual((self.again["ok"], self.again["did"]), (True, "nothing"),
                         self.again.get("error"))


@unittest.skipUnless(HAS_BACKENDS, "OpenMM, PDBFixer and MDTraj are needed")
class TestAStoppedExtensionDoesNotStopTheNext(unittest.TestCase):
    """The same, extended again by hand rather than carried on: the piece
    with no checkpoint is set aside and the new length runs from the study's
    own checkpoint, not refused for the piece left behind."""

    @classmethod
    def setUpClass(cls):
        from unittest import mock

        pytest.importorskip("openmm.app")
        import openmm.app

        from fastmdxplora import FastMDXplora
        from fastmdxplora.simulation.resume import extend_study, production_done_ns

        root = Path(tempfile.mkdtemp())
        cls.addClassCleanup(shutil.rmtree, root, True)
        pdb = root / "tri-ala.pdb"
        pdb.write_text(TRI_ALANINE)
        cls.study = root / "study"
        taken = {"steps": 0, "limit": 10 ** 9}
        with mock.patch.object(openmm.app.Simulation, "step", _stepping(taken)):
            FastMDXplora(config_data=_config(pdb, telemetry_interval=25),
                         output_dir=str(cls.study)).explore()
            taken.update(steps=0, limit=40)
            extend_study(cls.study, more_ns=0.0008)
        cls.answer = extend_study(cls.study, more_ns=0.0004)
        cls.done = production_done_ns(cls.study)

    def test_the_new_length_runs(self):
        self.assertTrue(self.answer["ok"], self.answer.get("error"))
        aside = self.study / "segment-001-stopped-before-its-checkpoint"
        self.assertEqual(self.answer.get("set_aside"), str(aside))
        self.assertAlmostEqual(self.done, 0.0010, places=6)


@unittest.skipUnless(HAS_BACKENDS, "OpenMM, PDBFixer and MDTraj are needed")
class TestAnExtensionStoppedPastItsCheckpointKeepsItsLength(unittest.TestCase):
    """Killed in production (its checkpoint at 200 steps), extended past its
    plan by 0.0008 ns and killed before the extension's checkpoint, carried
    on and killed again past the carry-on's checkpoint (250 of its 400
    steps): nothing of the plan is left, and the join was refused for the
    piece that did not finish at every resume after. The rest of the
    0.0008 ns asked is run."""

    @classmethod
    def setUpClass(cls):
        from unittest import mock

        pytest.importorskip("openmm.app")
        import openmm.app

        from fastmdxplora import FastMDXplora
        from fastmdxplora.simulation.resume import (
            extend_study,
            production_done_ns,
            resume_study,
        )

        root = Path(tempfile.mkdtemp())
        cls.addClassCleanup(shutil.rmtree, root, True)
        pdb = root / "tri-ala.pdb"
        pdb.write_text(TRI_ALANINE)
        cls.study = root / "study"
        taken = {"steps": 0, "limit": 200 + 250}
        with mock.patch.object(openmm.app.Simulation, "step", _stepping(taken)):
            FastMDXplora(config_data=_config(pdb, telemetry_interval=25),
                         output_dir=str(cls.study)).explore()
            (cls.study / "manifest.json").unlink()
            taken.update(steps=0, limit=40)
            extend_study(cls.study, more_ns=0.0008)
            taken.update(steps=0, limit=250)
            cls.stopped = resume_study(cls.study)
        cls.had_checkpoint = (cls.study / "segment-001" / "simulation" / "checkpoint.chk").is_file()
        cls.answer = resume_study(cls.study)
        cls.done = production_done_ns(cls.study)

    def test_the_carry_on_stopped_past_its_checkpoint(self):
        self.assertFalse(self.stopped["ok"])
        self.assertTrue(self.had_checkpoint)

    def test_the_rest_of_the_length_asked_runs(self):
        self.assertEqual((self.answer["ok"], self.answer["did"]), (True, "continued"),
                         self.answer.get("error"))
        # 200 steps to the first checkpoint, and the 0.0008 ns asked.
        self.assertAlmostEqual(self.done, 0.0012, places=6)


def _aside(root, name, *, asked, began, ran_in):
    folder = root / name
    folder.mkdir(parents=True)
    (folder / "resolved_config.yml").write_text(
        f"output: {ran_in}/segment-002\nsimulation:\n  duration_ns: {asked}\n"
        f"  resume_from: {began}\n", encoding="utf-8")
    return folder


class TestTheLengthAskedIsTheNewestSetAsideFromThisCheckpoint(unittest.TestCase):
    """Which set-aside piece's length a carry-on runs: one from the same
    checkpoint within the study (another piece's, or another study's, with
    the same name of file is not it), and of two from it, the one set aside
    last by its name, whatever its folder's time says after a copy."""

    def setUp(self):
        self.root = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, self.root, True)
        self.own = self.root / "simulation" / "checkpoint.chk"
        self.own.parent.mkdir(parents=True)
        self.own.write_bytes(b"")

    def test_another_piece_s_checkpoint_is_not_the_study_s(self):
        from fastmdxplora.simulation.resume import NEVER_CHECKPOINTED, _length_asked

        was = "/somewhere/old"
        _aside(self.root, f"segment-002{NEVER_CHECKPOINTED}", asked=50,
               began=f"{was}/segment-001/simulation/checkpoint.chk", ran_in=was)
        _aside(self.root, f"segment-003{NEVER_CHECKPOINTED}", asked=60,
               began="/another/study/simulation/checkpoint.chk", ran_in=was)
        self.assertIsNone(_length_asked(self.root, str(self.own)))
        _aside(self.root, f"segment-001{NEVER_CHECKPOINTED}", asked=7,
               began=f"{was}/simulation/checkpoint.chk", ran_in=was)
        self.assertEqual(_length_asked(self.root, str(self.own)), 7.0)

    def test_the_last_set_aside_by_its_name_is_read(self):
        import os

        from fastmdxplora.simulation.resume import NEVER_CHECKPOINTED, _length_asked

        first = _aside(self.root, f"segment-001{NEVER_CHECKPOINTED}", asked=20,
                       began=str(self.own), ran_in=str(self.root))
        second = _aside(self.root, f"segment-001{NEVER_CHECKPOINTED}-2", asked=5,
                        began=str(self.own), ran_in=str(self.root))
        # Copied without its times: the first now looks the newer.
        os.utime(second, (1_000_000, 1_000_000))
        os.utime(first, (2_000_000, 2_000_000))
        self.assertEqual(_length_asked(self.root, str(self.own)), 5.0)


class TestAPieceStillRunningIsNotSetAside(unittest.TestCase):
    """A resume while the study's extension is still going, before its
    checkpoint: refused as still running, its folder left where it is."""

    def test_it_is_refused_and_the_segment_left(self):
        from unittest import mock

        from fastmdxplora.simulation import resume

        root = _folder()
        self.addCleanup(shutil.rmtree, root, True)
        (root / "resolved_config.yml").write_text("simulation: {}\n", encoding="utf-8")
        segment = root / "segment-001"
        (segment / "simulation").mkdir(parents=True)
        with mock.patch.object(resume, "_still_running", lambda where: where == segment):
            answer = resume.resume_study(root)
        self.assertEqual((answer["ok"], answer["did"]), (False, "nothing"))
        self.assertIn("still going", answer["error"])
        self.assertTrue(segment.is_dir())
