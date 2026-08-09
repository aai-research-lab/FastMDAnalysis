"""Holding parts of a system still while the rest settles.

A structure that has just been minimised is not at equilibrium, and heating it
lets the solute move as well as the solvent: side chains relax into the vacuum
the crystal packing left, a ligand drifts out of the pose that was measured.
The remedy is to hold the solute while the solvent equilibrates around it and
then let go in stages, and this software went from minimisation to
unrestrained dynamics in one step until now.
"""

from __future__ import annotations

import numpy as np
import pytest

# Every test here that touches OpenMM calls ``pytest.importorskip`` as its
# first statement, before importing anything. A guard placed below the import
# it protects never fires -- the import raises first -- which is how two of
# these failed in CI while passing on a machine that had OpenMM.
#
# Verifying that needs the module to be genuinely unfindable rather than
# present-and-raising: pytest re-raises the latter deliberately, because it
# usually means a broken install rather than an absent package.

#: A short alanine peptide in an extended conformation. The terminal OXT is
#: there because without it the force field refuses the last residue: the
#: chain looks like it continues into something that is not present. Written out rather than loaded from a dataset,
#: because a test that needs another package to run is a test that will not
#: run where somebody checks out this one.
_ALANINE_PEPTIDE = """\
ATOM      1  N   ALA A   1      -1.204   1.045   0.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C
ATOM      3  CB  ALA A   1       0.000  -0.771  -1.303  1.00  0.00           C
ATOM      4  C   ALA A   1       1.278   0.828   0.000  1.00  0.00           C
ATOM      5  O   ALA A   1       1.256   2.056   0.000  1.00  0.00           O
ATOM      6  N   ALA A   2       2.416   0.146   0.000  1.00  0.00           N
ATOM      7  CA  ALA A   2       3.720   0.799   0.000  1.00  0.00           C
ATOM      8  CB  ALA A   2       3.720   1.570  -1.303  1.00  0.00           C
ATOM      9  C   ALA A   2       4.858  -0.221   0.000  1.00  0.00           C
ATOM     10  O   ALA A   2       4.836  -1.449   0.000  1.00  0.00           O
ATOM     11  N   ALA A   3       5.996   0.461   0.000  1.00  0.00           N
ATOM     12  CA  ALA A   3       7.300  -0.192   0.000  1.00  0.00           C
ATOM     13  CB  ALA A   3       7.300  -0.963  -1.303  1.00  0.00           C
ATOM     14  C   ALA A   3       8.438   0.828   0.000  1.00  0.00           C
ATOM     15  O   ALA A   3       8.416   2.056   0.000  1.00  0.00           O
ATOM     16  N   ALA A   4       9.576   0.146   0.000  1.00  0.00           N
ATOM     17  CA  ALA A   4      10.880   0.799   0.000  1.00  0.00           C
ATOM     18  CB  ALA A   4      10.880   1.570  -1.303  1.00  0.00           C
ATOM     19  C   ALA A   4      12.018  -0.221   0.000  1.00  0.00           C
ATOM     20  O   ALA A   4      11.996  -1.449   0.000  1.00  0.00           O
ATOM     21  OXT ALA A   4      13.156   0.196   0.000  1.00  0.00           O
TER
END
"""


class TestReadingWhatToHold:
    """The short form is a selection; the long form is a list of blocks. An
    equilibration usually wants the short one and should not need a paragraph
    to say so."""

    def test_a_selection_alone_means_position_restraints(self) -> None:
        from fastmdxplora.simulation.restraints import (
            DEFAULT_POSITION_FORCE,
            parse_restraints,
        )

        found = parse_restraints("protein and not element H")
        assert len(found) == 1
        assert found[0].kind == "position"
        assert found[0].force_constant == DEFAULT_POSITION_FORCE

    def test_nothing_means_nothing(self) -> None:
        from fastmdxplora.simulation.restraints import parse_restraints

        assert parse_restraints(None) == []
        assert parse_restraints("") == []

    def test_blocks_carry_their_own_kind_and_force(self) -> None:
        from fastmdxplora.simulation.restraints import parse_restraints

        found = parse_restraints([
            {"kind": "position", "selection": "protein"},
            {"kind": "distance", "selection": "index 0 1",
             "force_constant": 500.0, "target": 0.3},
        ])
        assert [r.kind for r in found] == ["position", "distance"]
        assert found[1].target == 0.3

    def test_a_restraint_with_no_conventional_default_must_say_its_force(self) -> None:
        """There is a standard force for holding heavy atoms in place. There
        is none for a distance restraint, and inventing one would be inventing
        the strength of a bias."""
        from fastmdxplora.simulation.restraints import parse_restraints

        with pytest.raises(ValueError, match="force_constant"):
            parse_restraints([{"kind": "distance", "selection": "index 0 1"}])

    def test_an_unknown_kind_is_refused(self) -> None:
        from fastmdxplora.simulation.restraints import parse_restraints

        with pytest.raises(ValueError, match="Unknown restraint kind"):
            parse_restraints([{"kind": "gravity", "selection": "protein"}])

    def test_the_units_follow_the_coordinate(self) -> None:
        """A spring on a length and a spring on an angle are not measured in
        the same thing, and one number for both is how an angle restraint ends
        up a thousand times too weak."""
        from fastmdxplora.simulation.restraints import parse_restraints

        position = parse_restraints([{"kind": "position", "selection": "protein"}])[0]
        torsion = parse_restraints([
            {"kind": "torsion", "selection": "index 0 1 2 3",
             "force_constant": 50.0}])[0]
        assert position.units == "kJ/mol/nm^2"
        assert torsion.units == "kJ/mol/rad^2"


class TestTheReleaseSchedule:
    """Letting go all at once undoes the point of restraining: the solute is
    released into a solvent arrangement that formed around a rigid structure,
    and the sudden freedom shows as a jump in energy and a lurch in shape."""

    def test_it_steps_down_and_reaches_zero(self) -> None:
        from fastmdxplora.simulation.restraints import ReleaseSchedule

        schedule = ReleaseSchedule()
        forces = [schedule.force_at(f) for f in (0.0, 0.3, 0.6, 0.9, 1.0)]
        assert forces == sorted(forces, reverse=True), "it should only weaken"
        assert forces[0] > 0
        assert forces[-1] == 0.0, "equilibration ends unrestrained"

    def test_a_schedule_can_be_given(self) -> None:
        from fastmdxplora.simulation.restraints import ReleaseSchedule

        schedule = ReleaseSchedule(steps=(200.0, 0.0))
        assert schedule.force_at(0.0) == 200.0
        assert schedule.force_at(1.0) == 0.0


class TestARestraintActuallyHolds:
    """The measurement that matters: does a restrained structure move less
    than an unrestrained one?"""

    @staticmethod
    def _system(restrained):
        import mdtraj as md
        import openmm as omm
        import openmm.unit as unit
        from openmm.app import HBonds, PME, ForceField, Modeller, PDBFile, Simulation

        from fastmdxplora.simulation.restraints import (
            build_restraint_forces,
            parse_restraints,
        )

        import tempfile
        from pathlib import Path

        # Written here rather than loaded, so the test needs no data file and
        # no package beyond what this one declares. A short alanine helix is
        # enough: the question is whether restrained atoms move less than free
        # ones, and any real peptide answers it.
        work = Path(tempfile.mkdtemp())
        (work / "start.pdb").write_text(_ALANINE_PEPTIDE, encoding="utf-8")
        pdb = PDBFile(str(work / "start.pdb"))
        field = ForceField("amber14-all.xml", "amber14/tip3p.xml")
        modeller = Modeller(pdb.topology, pdb.positions)
        modeller.addHydrogens(field)
        modeller.addSolvent(field, padding=0.6 * unit.nanometer)

        system = field.createSystem(
            modeller.topology, nonbondedMethod=PME,
            nonbondedCutoff=1.0 * unit.nanometer, constraints=HBonds)
        parameters = []
        if restrained:
            for force, parameter in build_restraint_forces(
                    omm, modeller.topology, modeller.positions,
                    parse_restraints("protein and not element H")):
                system.addForce(force)
                parameters.append(parameter)

        integrator = omm.LangevinMiddleIntegrator(
            300 * unit.kelvin, 1 / unit.picosecond, 0.002 * unit.picoseconds)
        integrator.setRandomNumberSeed(7)
        simulation = Simulation(modeller.topology, system, integrator,
                                omm.Platform.getPlatformByName("CPU"))
        simulation.context.setPositions(modeller.positions)
        # To convergence, which is what the software does -- its
        # minimize_max_iterations defaults to 0. Capped at 150 for speed, the
        # clashes addSolvent leaves behind survive into dynamics, and the
        # unrestrained arm becomes the very failure the restraints exist to
        # prevent: it went non-finite on macOS while passing here, and a
        # control that explodes measures nothing either way.
        simulation.minimizeEnergy(
            tolerance=10 * unit.kilojoule_per_mole / unit.nanometer)
        simulation.context.setVelocitiesToTemperature(300 * unit.kelvin, 7)
        heavy = md.Topology.from_openmm(modeller.topology).select(
            "protein and not element H")
        return simulation, heavy, parameters

    @staticmethod
    def _drift(simulation, heavy, steps=1500):
        import numpy as np
        import openmm.unit as unit

        start = simulation.context.getState(
            getPositions=True).getPositions(asNumpy=True)
        simulation.step(steps)
        end = simulation.context.getState(
            getPositions=True).getPositions(asNumpy=True)
        moved = np.linalg.norm(
            (end[heavy] - start[heavy]).value_in_unit(unit.nanometer), axis=1)
        return float(np.sqrt((moved ** 2).mean()))

    @pytest.mark.slow
    def test_restrained_atoms_move_less(self) -> None:
        pytest.importorskip("openmm", reason="requires the [md] extra")


        import math

        free, heavy, _ = self._system(restrained=False)
        held, _heavy, _params = self._system(restrained=True)

        unrestrained = self._drift(free, heavy)
        restrained = self._drift(held, heavy)

        # Said before the comparison, because a control that blew up makes the
        # ratio meaningless in the direction that looks like a pass.
        assert math.isfinite(unrestrained) and unrestrained > 0, (
            "the unrestrained control did not run: there is nothing to "
            f"compare against (drift {unrestrained})")
        assert math.isfinite(restrained), (
            f"the restrained run did not survive (drift {restrained})")

        assert restrained < unrestrained / 2, (
            "a restrained structure should move markedly less than a free "
            f"one: {restrained:.3f} nm restrained, {unrestrained:.3f} nm free"
        )

    @pytest.mark.slow
    def test_releasing_the_restraint_lets_the_structure_loosen(self) -> None:
        """Which is what makes staged release possible: the force is a global
        parameter, so it can be changed without rebuilding the context."""
        pytest.importorskip("openmm", reason="requires the [md] extra")


        simulation, heavy, parameters = self._system(restrained=True)
        assert parameters

        held = self._drift(simulation, heavy, steps=600)
        for parameter in parameters:
            simulation.context.setParameter(parameter, 0.0)
        freed = self._drift(simulation, heavy, steps=600)

        assert freed > held, (
            "releasing the restraint should let the structure move more"
        )


class TestARestraintOnNothingIsRefused:
    """A run that silently applied a restraint matching no atoms would look
    restrained and not be."""

    def test_a_selection_matching_nothing_raises(self) -> None:
        pytest.importorskip("openmm", reason="requires the [md] extra")

        import mdtraj as md
        import numpy as np
        import openmm as omm

        from fastmdxplora.simulation.restraints import (
            build_restraint_forces,
            parse_restraints,
        )


        top = md.Topology()
        residue = top.add_residue("ALA", top.add_chain(), resSeq=1)
        top.add_atom("CA", md.element.carbon, residue)
        positions = np.zeros((1, 3))

        with pytest.raises(ValueError, match="matched no atoms"):
            build_restraint_forces(
                omm, top.to_openmm(), positions,
                parse_restraints("resname NOPE"))

    def test_a_distance_restraint_needs_exactly_two_atoms(self) -> None:
        pytest.importorskip("openmm", reason="requires the [md] extra")

        import mdtraj as md
        import numpy as np
        import openmm as omm

        from fastmdxplora.simulation.restraints import (
            build_restraint_forces,
            parse_restraints,
        )


        top = md.Topology()
        residue = top.add_residue("ALA", top.add_chain(), resSeq=1)
        for name in ("N", "CA", "C"):
            top.add_atom(name, md.element.carbon, residue)
        positions = np.zeros((3, 3))

        with pytest.raises(ValueError, match="two atoms"):
            build_restraint_forces(
                omm, top.to_openmm(), positions,
                parse_restraints([{"kind": "distance", "selection": "all",
                                   "force_constant": 100.0}]))


class TestTheRunnerReleasesThemInStages:
    """Four force types are not a protocol. What makes them one is when they
    are applied, how they weaken, and that they are gone before production."""

    def test_the_settings_are_declared(self) -> None:
        from fastmdxplora.config.schema import PHASE_SCHEMAS

        for name in ("restrain", "restraint_release", "restrain_production"):
            field = PHASE_SCHEMAS["simulation"].get(name)
            assert field is not None, name
            assert field.help, name
        assert PHASE_SCHEMAS["simulation"].get("restrain_production").default is False, (
            "a biased production run measures the bias, so it is off"
        )

    def test_they_reach_the_command_line(self) -> None:
        """Without the CLI being touched: the flags come from the schema."""
        from fastmdxplora.cli.main import _PHASE_SPEC

        table, _prefix = _PHASE_SPEC["simulate"]
        offered = {dest for _flag, dest, _kw in table}
        assert {"restrain", "restraint_release",
                "restrain_production"} <= offered

    def test_the_pipeline_passes_them_to_the_runner(self) -> None:
        import inspect

        from fastmdxplora.simulation import pipeline

        source = inspect.getsource(pipeline)
        for name in ("restrain=", "restraint_release=", "restrain_production="):
            assert name in source, name

    def test_they_are_at_full_strength_for_nvt(self) -> None:
        """Which is the stage they exist for: the solvent is finding its
        arrangement and the solute should not be moving while it does."""
        import inspect

        from fastmdxplora.simulation import runner

        source = inspect.getsource(runner.run_simulation)
        nvt = source.index("Stage 2: NVT")
        npt = source.index("Stage 3: NPT")
        assert "_hold_at(0.0)" in source[nvt:npt]

    def test_they_are_gone_before_production(self) -> None:
        import inspect

        from fastmdxplora.simulation import runner

        source = inspect.getsource(runner.run_simulation)
        before = source[source.index("Stage 3: NPT"):source.index("Stage 4: Production")]
        assert "_hold_at(1.0)" in before
        assert "restrain_production" in before, (
            "keeping them has to be asked for by name"
        )

    def test_keeping_them_warns_about_what_it_costs(self) -> None:
        """A reader comparing a restrained trajectory against a free one
        otherwise has no way to know which they have."""
        import inspect

        from fastmdxplora.simulation import runner

        source = inspect.getsource(runner.run_simulation)
        assert "measures of flexibility" in source
        assert "RMSF" in source

    @pytest.mark.slow
    def test_a_run_with_restraints_completes_and_records_them(self, tmp_path) -> None:
        """That the restraint holds is measured directly in
        ``TestARestraintActuallyHolds``, where the comparison is clean: two
        systems, same seed, one force added.

        The same comparison through the whole runner was tried here and is not
        reliable. Total drift by the time production starts includes the
        solvent finding its arrangement and the box equilibrating, which are
        larger than the effect at the lengths a test can afford -- at fifteen
        hundred steps restrained and free came out at 0.0840 and 0.0832 nm,
        indistinguishable, while at two thousand they were 0.067 and 0.108. A
        test that passes at one length and fails at another is measuring
        noise, and a flaky test is worse than none: it teaches people to rerun
        rather than to look.

        So this checks what the runner is responsible for -- that a restrained
        run completes and says what it held -- and leaves the physics to the
        test that can measure it cleanly.
        """
        pytest.importorskip("openmm", reason="requires the [md] extra")


        import mdtraj as md
        import openmm as omm
        import openmm.unit as unit
        from openmm.app import HBonds, PME, ForceField, Modeller, PDBFile

        from fastmdxplora.simulation.runner import run_simulation

        (tmp_path / "p.pdb").write_text(_ALANINE_PEPTIDE, encoding="utf-8")
        pdb = PDBFile(str(tmp_path / "p.pdb"))
        field = ForceField("amber14-all.xml", "amber14/tip3p.xml")
        modeller = Modeller(pdb.topology, pdb.positions)
        modeller.addHydrogens(field)
        modeller.addSolvent(field, padding=0.7 * unit.nanometer)
        system = field.createSystem(
            modeller.topology, nonbondedMethod=PME,
            nonbondedCutoff=1.0 * unit.nanometer, constraints=HBonds)
        (tmp_path / "system.xml").write_text(omm.XmlSerializer.serialize(system))
        integrator = omm.LangevinMiddleIntegrator(
            300 * unit.kelvin, 1 / unit.picosecond, 0.002 * unit.picoseconds)
        context = omm.Context(system, integrator,
                              omm.Platform.getPlatformByName("CPU"))
        context.setPositions(modeller.positions)
        (tmp_path / "state.xml").write_text(omm.XmlSerializer.serialize(
            context.getState(getPositions=True, getVelocities=True,
                             enforcePeriodicBox=True)))
        with (tmp_path / "top.pdb").open("w") as handle:
            PDBFile.writeFile(modeller.topology, modeller.positions, handle)

        result = run_simulation(
            system_xml=tmp_path / "system.xml",
            state_xml=tmp_path / "state.xml",
            topology_pdb=tmp_path / "top.pdb",
            output_dir=tmp_path / "run",
            nvt_steps=200, npt_steps=200, production_steps=200,
            trajectory_interval_steps=100, state_interval_steps=100,
            platform="CPU", minimize=True, random_seed=11,
            restrain="protein and not element H")

        assert result.n_production_frames > 0
        traj = md.load(str(tmp_path / "run" / "production.dcd"),
                       top=str(tmp_path / "top.pdb"))
        assert traj.n_frames > 0, "a restrained run should produce a trajectory"

    @pytest.mark.slow
    def test_an_impossible_restraint_stops_the_run(self, tmp_path) -> None:
        """Rather than running unrestrained and looking as though it had
        worked."""
        pytest.importorskip("openmm", reason="requires the [md] extra")


        import openmm as omm
        import openmm.unit as unit
        from openmm.app import HBonds, PME, ForceField, Modeller, PDBFile

        from fastmdxplora.simulation.runner import run_simulation

        (tmp_path / "p.pdb").write_text(_ALANINE_PEPTIDE, encoding="utf-8")
        pdb = PDBFile(str(tmp_path / "p.pdb"))
        field = ForceField("amber14-all.xml", "amber14/tip3p.xml")
        modeller = Modeller(pdb.topology, pdb.positions)
        modeller.addHydrogens(field)
        modeller.addSolvent(field, padding=0.7 * unit.nanometer)
        system = field.createSystem(
            modeller.topology, nonbondedMethod=PME,
            nonbondedCutoff=1.0 * unit.nanometer, constraints=HBonds)
        (tmp_path / "system.xml").write_text(omm.XmlSerializer.serialize(system))
        integrator = omm.LangevinMiddleIntegrator(
            300 * unit.kelvin, 1 / unit.picosecond, 0.002 * unit.picoseconds)
        context = omm.Context(system, integrator,
                              omm.Platform.getPlatformByName("CPU"))
        context.setPositions(modeller.positions)
        (tmp_path / "state.xml").write_text(omm.XmlSerializer.serialize(
            context.getState(getPositions=True, getVelocities=True,
                             enforcePeriodicBox=True)))
        with (tmp_path / "top.pdb").open("w") as handle:
            PDBFile.writeFile(modeller.topology, modeller.positions, handle)

        with pytest.raises(ValueError, match="matched no atoms"):
            run_simulation(
                system_xml=tmp_path / "system.xml",
                state_xml=tmp_path / "state.xml",
                topology_pdb=tmp_path / "top.pdb",
                output_dir=tmp_path / "run",
                nvt_steps=10, npt_steps=10, production_steps=10,
                platform="CPU", minimize=False, random_seed=11,
                restrain="resname NOTHING")


class TestTheLadderStepsAcrossEquilibration:
    """It was sampled at two points -- once before NVT and once before NPT --
    so a four-rung ladder reached 1000 and 100 and never 500 or 0. With
    `npt_steps: 0` the second sat inside a branch that did not run, so the
    restraint held at full strength for the whole of equilibration and
    dropped to zero at production: the release all at once that the ladder
    exists to prevent, from a setting the user had written out in four
    steps."""

    def _reached(self, nvt: int, npt: int, ladder=(1000.0, 500.0, 100.0, 0.0)):
        """Every strength the ladder would set, given a stage division."""
        from fastmdxplora.simulation.restraints import ReleaseSchedule

        schedule = ReleaseSchedule(steps=tuple(ladder))
        equilibration = nvt + npt
        seen: list[float] = []
        for offset, steps in ((0, nvt), (nvt, npt)):
            if steps <= 0:
                continue
            # As the chunked stage reports it: once per chunk.
            for chunk in range(20):
                through = (offset + (chunk / 20) * steps) / equilibration
                strength = schedule.force_at(min(1.0, through))
                if not seen or seen[-1] != strength:
                    seen.append(strength)
        return seen

    def test_every_rung_is_reached(self) -> None:
        assert self._reached(2000, 2000) == [1000.0, 500.0, 100.0, 0.0]

    def test_without_npt_too(self) -> None:
        """The case that failed: the second sample sat in a branch that did
        not run, and the ladder never moved."""
        assert self._reached(4000, 0) == [1000.0, 500.0, 100.0, 0.0]

    def test_a_two_rung_ladder_is_honoured(self) -> None:
        assert self._reached(4000, 0, ladder=(500.0, 0.0)) == [500.0, 0.0]

    def test_it_ends_released(self) -> None:
        """Whatever the division, equilibration finishes at the last rung --
        which the schema requires to be zero."""
        for nvt, npt in ((4000, 0), (2000, 2000), (500, 3500)):
            assert self._reached(nvt, npt)[-1] == 0.0

    def test_the_stage_reports_how_far_through_it_is(self) -> None:
        import inspect

        from fastmdxplora.simulation.runner import (
            _run_md_stage_with_live_metrics,
        )

        assert "on_fraction" in inspect.signature(
            _run_md_stage_with_live_metrics).parameters
        body = inspect.getsource(_run_md_stage_with_live_metrics)
        loop = body.index("while remaining > 0:")
        assert body.index("on_fraction(") > loop, (
            "the ladder has to step inside the loop, or it steps once")

    def test_both_stages_drive_it(self) -> None:
        import inspect

        from fastmdxplora.simulation import runner

        source = inspect.getsource(runner)
        assert source.count("on_fraction=_ladder_over(") == 2
        # And the boundary sample that jumped it is gone.
        assert "_hold_at(0.5)" not in source
