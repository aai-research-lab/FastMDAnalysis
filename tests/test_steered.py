"""Pulling a system along a coordinate, and what that does and does not give.

Some things do not happen on their own within reach of a simulation. Steered
MD attaches a spring to a collective variable and moves the anchor, dragging
the system whether or not it wants to go.
"""

from __future__ import annotations

import pytest


def _topology():
    import mdtraj as md

    top = md.Topology()
    chain = top.add_chain()
    for index in range(30):
        residue = top.add_residue("ALA", chain, resSeq=index + 1)
        for name in ("N", "CA", "C", "O"):
            top.add_atom(name, md.element.carbon, residue)
    ligand_chain = top.add_chain()
    ligand = top.add_residue("BNZ", ligand_chain, resSeq=900)
    for index in range(6):
        top.add_atom(f"C{index}", md.element.carbon, ligand)
    return top


def _plan(**extra):
    from fastmdxplora.simulation.steered import plan_steered

    return plan_steered(dict({
        "collective_variable": "ligand_distance",
        "ligand_resname": "BNZ",
        "site_selection": "resid 1 to 3 and name CA",
        "to": 3.0,
    }, **extra), _topology())


class TestWhatItClaims:
    """A pathway and the work along it, not a free energy: the work depends on
    how fast the anchor moved, and a single fast pull overestimates a
    barrier.
    """

    def test_it_does_not_claim_a_free_energy(self) -> None:
        record = _plan(**{"from": 0.4}).as_record()
        assert "not a free energy" in record["gives"]
        assert "overestimates" in record["gives"]

    def test_the_module_says_what_it_is_for(self) -> None:
        """Generating starting structures for umbrella sampling, which does
        give a free energy from equilibrium sampling rather than from work
        done in a hurry."""
        import inspect

        from fastmdxplora.simulation import steered

        text = inspect.getdoc(steered).lower()
        assert "umbrella" in text
        assert "jarzynski" in text

    def test_it_reports_the_pulling_rate(self) -> None:
        """The number that decides whether the work means anything."""
        plan = _plan(**{"from": 0.4, "steps": 500000})
        assert plan.rate_per_ns(2.0) == pytest.approx(2.6, rel=0.01)

    def test_a_rate_needs_a_starting_value(self) -> None:
        """Without one the anchor starts wherever the system is, and how far
        it travels is not known in advance."""
        assert _plan(steps=500000).rate_per_ns(2.0) is None


class TestTheScriptItWrites:
    def test_it_is_a_moving_restraint(self) -> None:
        from fastmdxplora.simulation.steered import build_steered_script

        script = build_steered_script(_plan(**{"from": 0.4, "steps": 100000}))
        assert "MOVINGRESTRAINT" in script
        assert "AT0=0.4" in script and "AT1=3" in script
        assert "STEP1=100000" in script

    def test_it_records_the_work(self) -> None:
        """Without it a steered run has produced a trajectory and no number."""
        from fastmdxplora.simulation.steered import build_steered_script

        assert "pull.work" in build_steered_script(_plan())

    def test_the_script_says_the_rate_matters(self) -> None:
        from fastmdxplora.simulation.steered import build_steered_script

        script = build_steered_script(_plan())
        assert "overestimates a barrier" in script

    def test_plumed_counts_atoms_from_one(self) -> None:
        from fastmdxplora.simulation.steered import build_steered_script

        # The ligand is atoms 120-125 counting from zero.
        assert "ATOMS=121,122,123,124,125,126" in build_steered_script(_plan())


class TestWhatItRefuses:
    def test_a_pull_needs_a_destination(self) -> None:
        from fastmdxplora.simulation.steered import plan_steered

        with pytest.raises(ValueError, match="needs a `to`"):
            plan_steered({"collective_variable": "torsion",
                          "selection": "index 0 1 2 3"}, _topology())

    def test_steps_must_be_positive(self) -> None:
        with pytest.raises(ValueError, match="positive number of `steps`"):
            _plan(steps=0)

    def test_it_reuses_the_variables_metadynamics_offers(self) -> None:
        """Same five, resolved the same way, with the same refusals."""
        from fastmdxplora.simulation.metadynamics import COLLECTIVE_VARIABLES
        from fastmdxplora.simulation.steered import plan_steered

        with pytest.raises(ValueError, match="Unknown collective variable"):
            plan_steered({"collective_variable": "vibes", "to": 1.0},
                         _topology())
        assert len(COLLECTIVE_VARIABLES) == 5

    def test_an_unbounded_pull_is_not_refused(self) -> None:
        """Metadynamics refuses an unbounded ligand run because the bias fills
        a basin that never fills. A pull has a destination, so it ends."""
        assert _plan().to_value == 3.0


class TestItReachesTheRunner:
    def test_the_setting_is_declared(self) -> None:
        from fastmdxplora.config.schema import PHASE_SCHEMAS

        field = PHASE_SCHEMAS["simulation"].get("steered")
        assert field is not None
        assert "not a free energy" in field.help

    def test_it_reaches_the_command_line(self) -> None:
        from fastmdxplora.cli.main import _PHASE_SPEC

        table, _prefix = _PHASE_SPEC["simulate"]
        assert "steered" in {dest for _flag, dest, _kw in table}

    def test_the_pipeline_passes_it(self) -> None:
        import inspect

        from fastmdxplora.simulation import pipeline

        assert "steered=params.get" in inspect.getsource(pipeline)

    def test_steering_and_metadynamics_together_are_refused(self) -> None:
        """Two ways of moving the same coordinate, whose forces would add."""
        import inspect

        from fastmdxplora.simulation import runner

        source = inspect.getsource(runner.run_simulation)
        assert "steered and metadynamics" in source
        assert "their\n            \"forces would add" in source or \
               "forces would add" in source

    def test_the_script_is_written_where_it_can_be_read(self) -> None:
        import inspect

        from fastmdxplora.simulation import runner

        assert "steered.plumed" in inspect.getsource(runner.run_simulation)
