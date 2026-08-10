"""Metadynamics without writing PLUMED input.

PLUMED can express almost any enhanced-sampling scheme, and the cost of that
is a language to learn before running the commonest one. Most metadynamics on
a protein-ligand system biases one of a handful of things, and those do not
need a language.
"""

from __future__ import annotations

import pytest


def _topology():
    import mdtraj as md

    top = md.Topology()
    chain = top.add_chain()
    for index in range(6):
        residue = top.add_residue("ALA", chain, resSeq=index + 1)
        for name in ("N", "CA", "C", "O"):
            top.add_atom(name, md.element.carbon, residue)
    ligand_chain = top.add_chain()
    ligand = top.add_residue("BNZ", ligand_chain, resSeq=900)
    for index in range(6):
        top.add_atom(f"C{index}", md.element.carbon, ligand)
    return top


class TestChoosingWhatToBias:
    """Choosing a collective variable is the decision the method turns on: if
    it does not distinguish the states that matter, the surface converges and
    describes something that is not the system."""

    def test_every_variable_says_what_it_does_not_separate(self) -> None:
        """A name alone does not carry that, and it is the failure mode."""
        from fastmdxplora.simulation.metadynamics import COLLECTIVE_VARIABLES

        assert COLLECTIVE_VARIABLES
        for name, description in COLLECTIVE_VARIABLES.items():
            assert "Does not" in description or "does not" in description, name

    def test_an_unknown_variable_is_refused(self) -> None:
        from fastmdxplora.simulation.metadynamics import plan_from_config

        with pytest.raises(ValueError, match="Unknown collective variable"):
            plan_from_config({"collective_variable": "vibes", "sigma": 0.1},
                             _topology())

    def test_the_refusal_points_at_the_general_route(self) -> None:
        """Anything more elaborate is still written by hand: this is a shorter
        path to the common case, not a replacement for the general one."""
        from fastmdxplora.simulation.metadynamics import plan_from_config

        with pytest.raises(ValueError, match="PLUMED input directly"):
            plan_from_config({"collective_variable": "vibes", "sigma": 0.1},
                             _topology())


class TestTheHillWidthHasToBeGiven:
    """It should be about the size of the fluctuations within a single state,
    and there is no default that is right for an arbitrary coordinate."""

    def test_sigma_is_required(self) -> None:
        from fastmdxplora.simulation.metadynamics import plan_from_config

        with pytest.raises(ValueError, match="sigma"):
            plan_from_config(
                {"collective_variable": "radius_of_gyration"}, _topology())

    def test_and_the_refusal_says_what_a_reasonable_one_looks_like(self) -> None:
        from fastmdxplora.simulation.metadynamics import plan_from_config

        with pytest.raises(ValueError, match="0.05 nm"):
            plan_from_config(
                {"collective_variable": "radius_of_gyration"}, _topology())


class TestTheScriptItWrites:
    @staticmethod
    def _plan(**spec):
        from fastmdxplora.simulation.metadynamics import plan_from_config

        return plan_from_config(spec, _topology(), temperature_K=310.0)

    def test_plumed_counts_atoms_from_one(self) -> None:
        """Everything else here counts from zero. Getting it wrong biases the
        atom next door and produces a plausible surface for the wrong
        coordinate."""
        from fastmdxplora.simulation.metadynamics import build_plumed_script

        # Bounded, because an unbounded ligand run is refused: it would push
        # the ligand into bulk solvent and never come back.
        plan = self._plan(collective_variable="ligand_distance",
                          ligand_resname="BNZ",
                          site_selection="resid 1 to 3 and name CA",
                          sigma=0.05, walls={"upper": 2.5})
        script = build_plumed_script(plan)
        # The ligand is atoms 24-29 counting from zero.
        assert "ATOMS=25,26,27,28,29,30" in script

    def test_it_is_well_tempered_by_default(self) -> None:
        """Plain metadynamics deposits at full height forever, so the bias
        never settles and the surface never converges."""
        from fastmdxplora.simulation.metadynamics import build_plumed_script

        plan = self._plan(collective_variable="radius_of_gyration",
                          selection="protein and name CA", sigma=0.05)
        assert plan.bias_factor > 1.0
        assert "BIASFACTOR=" in build_plumed_script(plan)

    def test_it_records_the_variable_and_the_bias(self) -> None:
        """Without them there is no way to tell afterwards whether the run
        converged, and a run whose convergence cannot be checked has not
        measured a free energy."""
        from fastmdxplora.simulation.metadynamics import build_plumed_script

        plan = self._plan(collective_variable="radius_of_gyration",
                          selection="protein and name CA", sigma=0.05)
        script = build_plumed_script(plan)
        assert "PRINT ARG=cv,metad.bias" in script
        assert "FILE=COLVAR" in script
        assert "FILE=HILLS" in script

    def test_the_script_says_what_is_being_biased(self) -> None:
        """Somebody checking a result should be able to read what was biased
        without reading the module that wrote it."""
        from fastmdxplora.simulation.metadynamics import build_plumed_script

        plan = self._plan(collective_variable="torsion",
                          selection="resid 0 and name N CA C O", sigma=0.35)
        script = build_plumed_script(plan)
        assert script.startswith("#")
        assert "Does not" in script or "does not" in script

    def test_a_torsion_needs_exactly_four_atoms(self) -> None:
        """`resid 0` has exactly four in this topology, which is why the first
        version of this test selected four and asserted it would be
        refused."""
        with pytest.raises(ValueError, match="four atoms"):
            self._plan(collective_variable="torsion",
                       selection="resid 0 1", sigma=0.35)

    def test_a_ligand_distance_needs_somewhere_to_measure_to(self) -> None:
        with pytest.raises(ValueError, match="site_selection"):
            self._plan(collective_variable="ligand_distance",
                       ligand_resname="BNZ", sigma=0.05)

    def test_a_selection_matching_nothing_is_refused(self) -> None:
        with pytest.raises(ValueError, match="matched no atoms"):
            self._plan(collective_variable="radius_of_gyration",
                       selection="resname NOPE", sigma=0.05)

    def test_a_ligand_rmsd_needs_a_reference(self) -> None:
        """It is measured against a structure, and there is no sensible
        default for which one."""
        from fastmdxplora.simulation.metadynamics import build_plumed_script

        plan = self._plan(collective_variable="ligand_rmsd",
                          ligand_resname="BNZ", sigma=0.05,
                          walls={"upper": 1.0})
        with pytest.raises(ValueError, match="reference"):
            build_plumed_script(plan)


class TestWhatIsRecordedWithTheResults:
    def test_the_plan_carries_what_it_separates(self) -> None:
        from fastmdxplora.simulation.metadynamics import plan_from_config

        plan = plan_from_config(
            {"collective_variable": "radius_of_gyration",
             "selection": "protein and name CA", "sigma": 0.05},
            _topology())
        record = plan.as_record()
        assert record["well_tempered"] is True
        assert "does not separate" in record["what_it_separates"].lower() or \
               "Does not" in record["what_it_separates"]


class TestItReachesTheRunner:
    def test_the_setting_is_declared(self) -> None:
        from fastmdxplora.config.schema import PHASE_SCHEMAS

        field = PHASE_SCHEMAS["simulation"].get("metadynamics")
        assert field is not None
        assert "collective_variable" in field.help
        assert "sigma" in field.help

    def test_it_reaches_the_command_line(self) -> None:
        from fastmdxplora.cli.main import _PHASE_SPEC

        table, _prefix = _PHASE_SPEC["simulate"]
        assert "metadynamics" in {dest for _flag, dest, _kw in table}

    def test_the_pipeline_passes_it(self) -> None:
        import inspect

        from fastmdxplora.simulation import pipeline

        assert "metadynamics=params.get" in inspect.getsource(pipeline)

    def test_the_runner_writes_the_script_where_it_can_be_read(self) -> None:
        """It decides what the run measures, so it belongs beside the results
        rather than in memory."""
        import inspect

        from fastmdxplora.simulation import runner

        source = inspect.getsource(runner.run_simulation)
        assert "metadynamics.plumed" in source
        assert "script_path.write_text" in source

    def test_it_becomes_a_plumed_run(self) -> None:
        """The existing integration does the biasing; this is a shorter way to
        describe it, not a second mechanism."""
        import inspect

        from fastmdxplora.simulation import runner

        source = inspect.getsource(runner.run_simulation)
        block = source[source.index("if metadynamics:"):]
        assert 'plumed = {"enabled": True' in block[:2000]


class TestBoundingWhereTheLigandGoes:
    """A metadynamics run on a ligand's distance will, given time, push the
    ligand into bulk solvent -- where the landscape is flat and unbounded, so
    the bias fills a basin that is effectively infinite and the run never
    comes back to the question.

    Without a bound, a ligand-distance run is not wrong so much as
    unfinishable.
    """

    @staticmethod
    def _spec(**extra):
        return dict({
            "collective_variable": "ligand_distance",
            "ligand_resname": "BNZ",
            "site_selection": "resid 1 to 3 and name CA",
            "sigma": 0.05,
        }, **extra)

    def test_an_unbounded_ligand_run_is_refused(self) -> None:
        from fastmdxplora.simulation.metadynamics import plan_from_config

        with pytest.raises(ValueError, match="bulk solvent"):
            plan_from_config(self._spec(), _topology())

    def test_but_can_be_asked_for_deliberately(self) -> None:
        from fastmdxplora.simulation.metadynamics import plan_from_config

        plan = plan_from_config(self._spec(unbounded=True), _topology())
        assert plan.walls is None and plan.funnel is None

    def test_a_wall_bounds_how_far(self) -> None:
        from fastmdxplora.simulation.metadynamics import (
            build_plumed_script,
            plan_from_config,
        )

        plan = plan_from_config(
            self._spec(walls={"upper": 2.5}), _topology())
        script = build_plumed_script(plan)
        assert "UPPER_WALLS ARG=cv AT=2.5" in script

    def test_a_wall_with_neither_end_is_refused(self) -> None:
        from fastmdxplora.simulation.metadynamics import plan_from_config

        with pytest.raises(ValueError, match="wall nowhere"):
            plan_from_config(self._spec(walls={"kappa": 100}), _topology())

    def test_a_funnel_bounds_where(self) -> None:
        """A flat wall bounds how far a ligand goes and not where, so the run
        still explores a whole shell of unbound positions at that distance."""
        from fastmdxplora.simulation.metadynamics import (
            build_plumed_script,
            plan_from_config,
        )

        plan = plan_from_config(
            self._spec(funnel={"axis_selection": "resid 4 5 and name CA"}),
            _topology())
        script = build_plumed_script(plan)
        assert "UPPER_WALLS ARG=outside AT=0" in script
        assert "proj: CUSTOM" in script and "rad: CUSTOM" in script

    def test_the_funnel_uses_only_standard_plumed(self) -> None:
        """Rather than the FUNNEL module, which is not in every build."""
        from fastmdxplora.simulation.metadynamics import (
            build_plumed_script,
            plan_from_config,
        )

        plan = plan_from_config(
            self._spec(funnel={"axis_selection": "resid 4 5 and name CA"}),
            _topology())
        script = build_plumed_script(plan)
        assert "FUNNEL " not in script.replace("# Funnel", "")
        for action in ("COM", "DISTANCE", "CUSTOM", "UPPER_WALLS"):
            assert action in script

    def test_it_is_wide_at_the_site_and_narrow_in_bulk(self) -> None:
        """Limongelli's shape: a cone over the site mouth opening into a
        cylinder, so the unbound state has a defined volume -- which is what
        makes an absolute binding free energy recoverable."""
        from fastmdxplora.simulation.metadynamics import (
            build_plumed_script,
            plan_from_config,
        )

        plan = plan_from_config(
            self._spec(funnel={"axis_selection": "resid 4 5 and name CA",
                               "cylinder_radius_nm": 0.1,
                               "switch_distance_nm": 1.5,
                               "alpha_rad": 0.55}),
            _topology())
        line = next(l for l in build_plumed_script(plan).splitlines()
                    if l.startswith("limit:"))
        assert "max(0.1," in line, "it never narrows below the cylinder"
        assert "(1.5-p)" in line, "and widens back towards the site"

    def test_a_funnel_needs_the_direction_the_ligand_leaves_by(self) -> None:
        """Nothing here can work that out, and one pointed the wrong way
        blocks the exit instead of following it."""
        from fastmdxplora.simulation.metadynamics import plan_from_config

        with pytest.raises(ValueError, match="axis_selection"):
            plan_from_config(self._spec(funnel={}), _topology())

    def test_a_funnel_only_applies_to_a_ligand_leaving(self) -> None:
        from fastmdxplora.simulation.metadynamics import plan_from_config

        with pytest.raises(ValueError, match="applies to ligand_distance"):
            plan_from_config(
                {"collective_variable": "radius_of_gyration",
                 "selection": "protein and name CA", "sigma": 0.05,
                 "funnel": {"axis_selection": "resid 4 and name CA"}},
                _topology())

    def test_a_bounded_run_records_its_bounds(self) -> None:
        from fastmdxplora.simulation.metadynamics import plan_from_config

        plan = plan_from_config(self._spec(walls={"upper": 2.5}), _topology())
        assert plan.as_record()["walls"]["upper"] == 2.5


class TestTheThreeAddedVariables:
    """Coordination, membrane depth and angle -- the ones missing from the
    first five that people actually use."""

    @staticmethod
    def _membrane_topology():
        import mdtraj as md

        top = _topology()
        lipids = top.add_chain()
        for index in range(8):
            lipid = top.add_residue("POP", lipids, resSeq=500 + index)
            for atom in range(4):
                top.add_atom(f"C{atom}", md.element.carbon, lipid)
        return top

    def test_coordination_counts_through_a_switching_function(self) -> None:
        """A hard cutoff has an infinite derivative at the boundary, and a
        bias needs a force."""
        from fastmdxplora.simulation.metadynamics import (
            build_plumed_script,
            plan_from_config,
        )

        plan = plan_from_config({
            "collective_variable": "coordination",
            "selection_a": "resname BNZ",
            "selection_b": "resid 1 to 3 and name CA",
            "sigma": 0.5}, _topology())
        script = build_plumed_script(plan)
        assert "COORDINATION" in script
        assert "NN=6 MM=12" in script, "a switching function, not a step"

    def test_the_switching_distance_is_a_setting_that_travels(self) -> None:
        """It was declared on the plan and read by the script and never
        passed from the config, so every run used the default whatever was
        asked for."""
        from fastmdxplora.simulation.metadynamics import (
            build_plumed_script,
            plan_from_config,
        )

        for wanted in (0.3, 0.45):
            plan = plan_from_config({
                "collective_variable": "coordination",
                "selection_a": "resname BNZ",
                "selection_b": "resid 1 to 3 and name CA",
                "sigma": 0.5, "coordination_r0": wanted}, _topology())
            assert plan.coordination_r0 == wanted
            assert f"R_0={wanted:g}" in build_plumed_script(plan)

    def test_coordination_needs_two_groups(self) -> None:
        from fastmdxplora.simulation.metadynamics import plan_from_config

        with pytest.raises(ValueError, match="selection_b"):
            plan_from_config({
                "collective_variable": "coordination",
                "selection_a": "resname BNZ", "sigma": 0.5}, _topology())

    def test_membrane_depth_is_measured_from_the_bilayer_itself(self) -> None:
        """A membrane drifts in the box over a long run, so depth against a
        fixed plane slowly becomes depth against nothing."""
        from fastmdxplora.simulation.metadynamics import (
            build_plumed_script,
            plan_from_config,
        )

        plan = plan_from_config({
            "collective_variable": "membrane_depth",
            "ligand_resname": "BNZ", "sigma": 0.05},
            self._membrane_topology())
        script = build_plumed_script(plan)
        assert "mem: COM" in script, "the bilayer's own centre"
        assert "sep.z" in script, "along the normal"

    def test_membrane_depth_needs_a_membrane(self) -> None:
        from fastmdxplora.simulation.metadynamics import plan_from_config

        with pytest.raises(ValueError, match="matched no atoms"):
            plan_from_config({
                "collective_variable": "membrane_depth",
                "ligand_resname": "BNZ", "sigma": 0.05}, _topology())

    def test_an_angle_is_three_atoms(self) -> None:
        from fastmdxplora.simulation.metadynamics import (
            build_plumed_script,
            plan_from_config,
        )

        plan = plan_from_config({
            "collective_variable": "angle",
            "selection": "resid 0 and name N CA C", "sigma": 0.2}, _topology())
        assert "ANGLE ATOMS=" in build_plumed_script(plan)

        with pytest.raises(ValueError, match="three atoms"):
            plan_from_config({
                "collective_variable": "angle",
                "selection": "resid 0", "sigma": 0.2}, _topology())

    def test_all_eight_say_what_they_do_not_separate(self) -> None:
        from fastmdxplora.simulation.metadynamics import COLLECTIVE_VARIABLES

        assert len(COLLECTIVE_VARIABLES) == 8
        for name, description in COLLECTIVE_VARIABLES.items():
            assert "oes not" in description, name

    def test_steered_md_gets_them_too(self) -> None:
        """It reuses the same layer, so a new variable arrives in both."""
        from fastmdxplora.simulation.steered import build_steered_script, plan_steered

        plan = plan_steered({
            "collective_variable": "membrane_depth",
            # Zero is the bilayer midplane here -- a real depth, not an
            # absent one -- which is why a pull has to say where it starts.
            "ligand_resname": "BNZ", "from": 0.0, "to": 2.5},
            self._membrane_topology())
        assert "sep.z" in build_steered_script(plan)


def test_one_translation_serves_every_method_that_biases() -> None:
    """Steered MD had a copy of the collective-variable translation, so
    membrane_depth reached metadynamics and fell through steering's else
    branch into one expecting a radius of gyration.

    A variable should arrive everywhere at once. Checked by counting where
    the translation lives rather than by listing the methods, because the
    list is the thing that drifts.
    """
    import inspect

    from fastmdxplora.simulation import metadynamics, steered, umbrella

    definitions = sum(
        inspect.getsource(module).count("cv: COORDINATION")
        for module in (metadynamics, steered, umbrella)
    )
    assert definitions == 1, (
        "the translation from a named variable to PLUMED should exist once"
    )

    # And steering reaches it rather than reimplementing.
    assert "cv_lines(" in inspect.getsource(steered.build_steered_script)


class TestTheStagePlanSurvivesTheBiasPlan:
    """`plan = plan_from_config(...)` in the runner rebound the stage plan --
    a dict of step counts that the rest of the function subscripts -- to a
    MetadynamicsPlan. The next `plan["nvt_steps"]` raised
    "'MetadynamicsPlan' object is not subscriptable", so metadynamics failed
    on the first line of real use.

    Thirty-eight tests cover this module and twenty-three the surface it
    produces. None of them runs the runner, which is where the two names met.
    """

    def _runner_source(self) -> str:
        import inspect

        from fastmdxplora.simulation import runner

        return inspect.getsource(runner)

    def test_the_bias_plan_has_its_own_name(self) -> None:
        source = self._runner_source()
        assert "bias_plan = plan_from_config(" in source
        assert "\n        plan = plan_from_config(" not in source

    def test_the_umbrella_path_already_did_this(self) -> None:
        """`cv_plan` beside it: the pattern was there to copy."""
        source = self._runner_source()
        assert "cv_plan" in source

    def test_the_stage_plan_is_still_a_mapping_afterwards(self) -> None:
        """The failure was a dict becoming a dataclass, so the guard is that
        every use of `plan` in the runner subscripts it."""
        import re

        source = self._runner_source()
        # Every `plan = ` assignment in the runner should be a mapping, so no
        # assignment from a *_from_config factory may bind the bare name.
        for match in re.finditer(r"^\s+plan = (\w+)\(", source, re.M):
            assert not match.group(1).endswith("_from_config"), match.group(0)
