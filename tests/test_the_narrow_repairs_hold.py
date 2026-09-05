"""Batch 4 of the 2026-09-04 audit: eight small things, each independent.

Nothing here is a large mechanism. What they share is that each one was a
decision taken correctly in one place and not carried to the next: the
clustering module's pairwise-RMSD helper takes an atom selection and the
dimensionality-reduction module's copy did not; `superposed()` fixed a whole
class of bug and `explore()` still accumulated its results; `describe_pmf`
was written to stop a barrier being read across ground no window visited and
`closure_gap` beside it read across ground no coordinate can occupy.
"""

from __future__ import annotations

import mdtraj as md
import numpy as np
import pytest


def _states_plus_solvent(seed: int = 0):
    """Twelve CA atoms in two conformational states, plus diffusing water.

    The contrast between states is large on the CA atoms and swamped by the
    solvent, which is what makes this the fixture the MDS bug needed and the
    solvent-free peptide in the existing test could not be.
    """
    rng = np.random.RandomState(seed)
    top = md.Topology()
    chain = top.add_chain()
    for _ in range(12):
        residue = top.add_residue("ALA", chain)
        top.add_atom("CA", md.element.carbon, residue)
    for _ in range(200):
        water = top.add_residue("HOH", chain)
        top.add_atom("O", md.element.oxygen, water)

    along = np.linspace(0.0, 1.1, 12)
    frames = []
    for frame in range(20):
        # A hinge, not a translation: `md.rmsd` superposes before measuring,
        # so a rigidly shifted copy of the same shape reads as identical and
        # would have made this fixture a fixture of nothing.
        bend = 0.0 if frame < 10 else 0.9
        state = np.stack(
            [along, bend * np.clip(along - 0.55, 0, None) ** 2 * 4.0,
             np.zeros_like(along)], axis=1)
        state = state + rng.normal(0, 0.01, state.shape)
        solvent = rng.uniform(0, 4.0, (200, 3))
        frames.append(np.vstack([state, solvent]))
    return md.Trajectory(np.array(frames, dtype=np.float32), top)


class TestDimensionalityReductionMeasuresTheSelection:
    """AUD13. `cluster._pairwise_rmsd` takes `atom_idx`; dimred's did not."""

    def test_the_helper_takes_a_selection(self) -> None:
        import inspect

        from fastmdxplora.analysis.dimred import _pairwise_rmsd

        assert "atom_indices" in inspect.signature(_pairwise_rmsd).parameters

    def test_the_selection_is_what_separates_the_states(self) -> None:
        from fastmdxplora.analysis.dimred import _pairwise_rmsd

        traj = _states_plus_solvent()
        alpha = traj.topology.select("name CA")

        selected = _pairwise_rmsd(traj, alpha)
        everything = _pairwise_rmsd(traj)

        def contrast(matrix):
            between = matrix[:10, 10:].mean()
            within = (matrix[:10, :10].mean() + matrix[10:, 10:].mean()) / 2
            return between / max(within, 1e-9)

        assert contrast(selected) > 5.0
        assert contrast(everything) < 2.0

    def test_mds_and_pca_no_longer_disagree(self) -> None:
        """The same call resolved the transition with PCA and hid it with
        MDS, which is how a difference of this size stays invisible: one of
        the two methods always looked right."""
        from fastmdxplora.analysis.dimred import DimRed

        traj = _states_plus_solvent()

        def separated(method):
            embedding = DimRed(
                methods=[method], selection="name CA").compute(traj)[method]
            first, second = embedding[:10], embedding[10:]
            gap = float(np.linalg.norm(first.mean(0) - second.mean(0)))
            spread = float(np.std(first) + np.std(second))
            return gap > spread

        assert separated("pca")
        assert separated("mds")


class TestEveryHistidineKeepsItsChemistry:
    """AUD9. An imidazole is aromatic whatever the force field calls it."""

    @pytest.mark.parametrize(
        "resname", ["HIS", "HIE", "HID", "HIP", "HSD", "HSE", "HSP"])
    def test_it_has_a_ring(self, resname: str) -> None:
        from fastmdxplora.analysis.interactions import _AROMATIC_RINGS

        assert resname in _AROMATIC_RINGS

    @pytest.mark.parametrize("resname", ["HIP", "HSP"])
    def test_the_doubly_protonated_ones_are_positive(self, resname) -> None:
        """The residue name *is* the setup phase's protonation decision, so
        leaving them out did not defer to it -- it discarded it."""
        from fastmdxplora.analysis.interactions import _POSITIVE_GROUPS

        assert resname in _POSITIVE_GROUPS

    @pytest.mark.parametrize("resname", ["HIS", "HIE", "HID", "HSD", "HSE"])
    def test_the_singly_protonated_ones_are_not(self, resname) -> None:
        """Still deferred to the setup phase, as the table's docstring says."""
        from fastmdxplora.analysis.interactions import _POSITIVE_GROUPS

        assert resname not in _POSITIVE_GROUPS

    @pytest.mark.parametrize("resname", ["HIE", "HID", "HIP"])
    def test_a_ring_is_found_on_a_real_topology(self, resname: str) -> None:
        from fastmdxplora.analysis.interactions import protein_aromatic_rings

        top = md.Topology()
        residue = top.add_residue(resname, top.add_chain())
        for name in ("CG", "ND1", "CD2", "CE1", "NE2"):
            element = md.element.nitrogen if name.startswith("N") else md.element.carbon
            top.add_atom(name, element, residue)

        rings = protein_aromatic_rings(top, range(top.n_atoms))
        assert len(rings) == 1 and len(rings[0]) == 5

    def test_they_are_no_longer_called_uninteresting(self) -> None:
        """They were on the list of residues with "nothing to contribute",
        so an unexamined histidine was not reported as unexamined either."""
        from fastmdxplora.analysis.interactions import residues_not_covered

        top = md.Topology()
        residue = top.add_residue("HIE", top.add_chain())
        for name in ("CG", "ND1", "CD2", "CE1", "NE2"):
            element = md.element.nitrogen if name.startswith("N") else md.element.carbon
            top.add_atom(name, element, residue)

        # Covered because it is in the ring table, not because it was
        # declared dull.
        assert residues_not_covered(top, range(top.n_atoms)) == {}


class TestASelectionThatDropsProteinSaysSo:
    """AUD9's other half: MDTraj's `protein` excludes HIE, HID and HSP."""

    def _peptide(self, middle: str):
        top = md.Topology()
        chain = top.add_chain()
        for name in ("ALA", middle, "ALA"):
            residue = top.add_residue(name, chain)
            for atom, element in (("N", "N"), ("CA", "C"),
                                  ("C", "C"), ("O", "O")):
                top.add_atom(atom, md.element.get_by_symbol(element), residue)
        xyz = np.zeros((2, top.n_atoms, 3), dtype=np.float32)
        xyz[:] = np.arange(top.n_atoms)[None, :, None] * 0.15
        return md.Trajectory(xyz, top)

    @pytest.mark.parametrize("resname", ["HIE", "HID"])
    def test_the_finding_names_the_residue(self, resname: str) -> None:
        from fastmdxplora.analysis.sasa import SASA

        analysis = SASA()
        analysis.select_atoms(self._peptide(resname))

        note = analysis.findings.get("selection_dropped_residues")
        assert note is not None, (
            f"{resname} is outside MDTraj's 'protein' and nothing said so"
        )
        assert resname in note

    def test_nothing_is_said_when_nothing_is_dropped(self) -> None:
        from fastmdxplora.analysis.sasa import SASA

        analysis = SASA()
        analysis.select_atoms(self._peptide("HIS"))
        assert "selection_dropped_residues" not in analysis.findings


class TestABondAngleIsNotACircle:
    """AUD20. PLUMED's ANGLE has domain [0, pi]; TORSION wraps."""

    def test_only_a_torsion_is_periodic(self) -> None:
        from fastmdxplora.simulation.umbrella import PERIODIC_VARIABLES

        assert "torsion" in PERIODIC_VARIABLES
        assert "angle" not in PERIODIC_VARIABLES

    def test_the_metadynamics_side_agrees(self) -> None:
        import inspect

        from fastmdxplora.simulation import metad_surface

        source = inspect.getsource(metad_surface.periodic_dimensions)
        circular = source[source.index("circular = ("):]
        circular = circular[:circular.index(")") + 1]
        assert "TORSION" in circular
        assert "ANGLE" not in circular


class TestAClosureGapNeedsSomethingToClose:
    """AUD21. `closure_gap` judged from the span alone."""

    def _profile(self, low, high):
        coordinate = np.linspace(low, high, 60)
        # A single well: the two ends differ, which is what closure_gap
        # measures and what makes a fabricated value visible.
        energy = 20.0 * (coordinate - low) / (high - low)
        return coordinate, energy

    def test_a_distance_pmf_gets_none(self) -> None:
        """0.3 to 6.8 nm spans more than 2*pi, which was the whole test."""
        from fastmdxplora.simulation.umbrella import describe_pmf

        coordinate, energy = self._profile(0.3, 6.8)
        summary = describe_pmf(coordinate, energy, periodic=False)

        assert summary["closure_gap_kjmol"] is None

    def test_a_torsion_still_gets_one(self) -> None:
        from fastmdxplora.simulation.umbrella import describe_pmf

        coordinate, energy = self._profile(-np.pi, np.pi)
        summary = describe_pmf(coordinate, energy, periodic=True)

        assert summary["closure_gap_kjmol"] is not None


class TestASettingThatDoesNothingIsRefused:
    """AUD22. `equilibration_steps` was accepted, recorded, and read by
    nothing -- the same failure as a typo, with better spelling."""

    def test_it_is_refused_by_name(self) -> None:
        from fastmdxplora.config.loader import ConfigError
        from fastmdxplora.simulation.umbrella import check_umbrella_keys

        with pytest.raises(ConfigError) as refusal:
            check_umbrella_keys({"equilibration_steps": 500_000})
        assert "equilibration_fraction" in str(refusal.value)

    def test_the_setting_that_works_is_still_accepted(self) -> None:
        from fastmdxplora.simulation.umbrella import check_umbrella_keys

        check_umbrella_keys({"equilibration_fraction": 0.3})

    def test_it_is_gone_from_the_record(self) -> None:
        """`pmf.json` used to carry it, which said a discard had happened."""
        from fastmdxplora.simulation.umbrella import plan_windows

        plan = plan_windows({
            "collective_variable": "distance",
            "selection_a": "resid 1", "selection_b": "resid 2",
            "from": 0.3, "to": 0.9, "n_windows": 3, "force_constant": 1000,
        })
        assert "equilibration_steps" not in plan.as_record()


class TestTheSecondCallIsNotToldAboutTheFirst:
    """AUD24. `explore()` appended to results that __init__ had set."""

    def test_a_retry_after_a_failure_reports_the_retry(self, tmp_path) -> None:
        from fastmdxplora import FastMDXplora
        from fastmdxplora.orchestrator import PhaseResult

        study = FastMDXplora(system="1UBQ", output_dir=tmp_path / "run")
        failing = {"n": 0}

        def _phase(phase, options):
            failing["n"] += 1
            status = "error" if failing["n"] == 1 else "ok"
            return PhaseResult(name=phase, status=status,
                               message="", output_dir=study.output_dir)

        study._run_phase = _phase
        study._write_manifest = lambda *a, **k: None

        first = study.explore(include=["setup"], report=False)[0]
        second = study.explore(include=["setup"], report=False)[0]

        assert first.status == "error"
        assert second.status == "ok", (
            "the second call inherited the first call's failure"
        )
        assert [p.status for p in second.phases] == ["ok"], (
            "and its phase list still carries the first call's records"
        )


class TestTheReportersAreClosedWhateverHappens:
    """AUD25. `_detach_all_reporters` was the last statement in the try."""

    def test_it_is_in_the_finally(self) -> None:
        """Read from the source, and stated as such.

        The behaviour needs a real OpenMM simulation raising mid-production,
        which this environment cannot always build. What can be checked
        cheaply is the structural fact that made the leak possible: the
        detach sitting where an exception skips it.
        """
        import inspect

        from fastmdxplora.simulation import runner

        source = inspect.getsource(runner.run_simulation)
        detach = source.index("_detach_all_reporters(simulation)")
        finally_at = source.rindex("finally:", 0, detach)
        except_at = source.rindex("except Exception as exc:", 0, detach)

        assert finally_at > except_at, (
            "the detach is inside the try or the except, so an exception "
            "between attaching a reporter and reaching it leaks the handle"
        )
