"""Every rule that measures a distance measures it under the minimum image.

The TODO records this class of defect three times -- `pl_interactions` and
`pl_contacts` through a stale unit cell (fixed in `a689c7e`), then
`ligand_rmsd` reading no cell at all. It is nine. Six interaction rules
ignored the box entirely: four passed a hardcoded `False` to the distance
helper and two took raw Cartesian separation between group centres, and
none of the six accepted a `periodic` parameter, so a user who knew could
not opt in -- while `options.json` recorded the treatment as applied.

They survived because **every one of the 32 tests in
`test_what_holds_the_ligand_is_measured.py` builds a trajectory with no
unit cell**, which makes `periodic` a no-op in all of them. A fixture that
cannot express the failure cannot catch it. So these fixtures carry a box,
and every partner in them sits across a face.

The geometry is deliberately trivial: two atoms 1.5 A apart through the
x face of a 5 nm cube read as 48.5 A when the box is ignored, and every
threshold in this module is between 2.5 and 6 A.
"""

from __future__ import annotations

import mdtraj as md
import numpy as np
import pytest

from fastmdxplora.analysis import interactions as rules

BOX_NM = 5.0
#: Just inside the face at each end, so the pair is 0.15 nm apart through
#: it and 4.85 nm apart across the middle of the box.
NEAR_LOW = 0.05
NEAR_HIGH = 4.90


def _topology(spec):
    """A topology from ``[(resname, [(name, element), ...]), ...]``."""
    top = md.Topology()
    chain = top.add_chain()
    for resname, atoms in spec:
        residue = top.add_residue(resname, chain)
        for name, element in atoms:
            top.add_atom(name, md.element.get_by_symbol(element), residue)
    return top


def _trajectory(top, xyz, *, box=BOX_NM):
    traj = md.Trajectory(np.array([xyz], dtype=np.float32), top)
    traj.unitcell_vectors = np.array(
        [np.eye(3) * box], dtype=np.float32)
    return traj


def _ring(centre_x, *, y=2.5, radius=0.14):
    """Six coplanar carbons, as a benzene ring lying in the xy plane."""
    angles = np.linspace(0.0, 2.0 * np.pi, 6, endpoint=False)
    return [[centre_x + radius * np.cos(a), y + radius * np.sin(a), 2.5]
            for a in angles]


def _chemistry(smiles: str, resname: str = "LIG"):
    """The same construction the sibling interaction suite uses."""
    from rdkit import Chem

    from fastmdxplora.analysis.ligand_chemistry import ResolvedChemistry

    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    return ResolvedChemistry(
        mol=mol, source="supplied", detail="built for a test",
        resname=resname, n_atoms=mol.GetNumAtoms(),
        charge_was_ambiguous=False)


class TestEveryRuleAcceptsThePeriodicOption:
    """The parameter must exist on all eight before it can be given to all
    eight. Six did not have it."""

    @pytest.mark.parametrize("name", [
        "hydrogen_bonds", "hydrophobic_contacts", "salt_bridges",
        "pi_stacking", "pi_cation", "halogen_bonds", "metal_coordination",
        "water_bridges",
    ])
    def test_it_takes_periodic(self, name: str) -> None:
        import inspect

        signature = inspect.signature(getattr(rules, name))
        assert "periodic" in signature.parameters, (
            f"{name} cannot be told about the unit cell, so the documented "
            f"`periodic` option cannot reach it."
        )

    @pytest.mark.parametrize("name", [
        "salt_bridges", "pi_stacking", "pi_cation", "halogen_bonds",
        "metal_coordination", "water_bridges",
    ])
    def test_the_default_is_to_use_the_box(self, name: str) -> None:
        """A trajectory that carries a cell is periodic; measuring it as if
        it were not is the unusual request, so it is the one that must be
        asked for."""
        import inspect

        assert inspect.signature(
            getattr(rules, name)).parameters["periodic"].default is True


class TestTheRuleIsGivenTheOption:

    def test_every_dispatch_line_passes_periodic(self) -> None:
        """Line-anchored rather than a substring over the whole method.

        This is a wiring check and its limits are worth stating: it reads
        source, which is the weaker kind of test, and it would pass on a
        dispatch that never ran. What makes it worth having is that the
        classes below measure the geometry through each rule directly, so
        this only has to answer the one question they cannot -- whether the
        *option* reaches all eight. It reached two.
        """
        import inspect

        from fastmdxplora.analysis import pl_interactions

        source = inspect.getsource(pl_interactions.ProteinLigandInteractions)
        lines = source.splitlines()
        starts = [i for i, line in enumerate(lines)
                  if line.strip().startswith('attempt("')]
        assert len(starts) == 8, f"expected 8 rules dispatched, found {len(starts)}"

        without = []
        for index in starts:
            # A dispatch is one call spread over two lines.
            call = " ".join(lines[index:index + 3])
            kind = call.split('"')[1]
            if "periodic=self.periodic" not in call:
                without.append(kind)

        assert not without, (
            f"dispatched without the box: {without}. `options.json` records "
            f"`periodic` as applied to all of them."
        )


class TestAContactThroughAFaceIsFound:
    """One geometry, eight rules, each with the partners across the box."""

    def test_salt_bridge(self) -> None:
        top = _topology([("LIG", [("N1", "N")]), ("ASP", [("OD1", "O"),
                                                          ("OD2", "O")])])
        traj = _trajectory(top, [
            [NEAR_HIGH, 2.5, 2.5],
            [NEAR_LOW, 2.5, 2.5],
            [NEAR_LOW, 2.5, 2.5],
        ])
        chemistry = _chemistry("[NH4+]")

        across = rules.salt_bridges(traj, chemistry, [0], [1, 2])
        ignored = rules.salt_bridges(
            traj, chemistry, [0], [1, 2], periodic=False)

        assert len(across) == 1
        assert across[0].distance_nm == pytest.approx(0.15, abs=0.02)
        assert ignored == []

    def test_metal_coordination(self) -> None:
        top = _topology([("ZN", [("ZN", "Zn")]), ("HIS", [("NE2", "N")])])
        traj = _trajectory(top, [[NEAR_HIGH, 2.5, 2.5], [NEAR_LOW, 2.5, 2.5]])

        across = rules.metal_coordination(traj, [0], [1])
        ignored = rules.metal_coordination(traj, [0], [1], periodic=False)

        assert len(across) == 1
        assert ignored == []

    def test_pi_stacking(self) -> None:
        """Two parallel rings 0.18 nm apart through the face."""
        top = _topology([
            ("LIG", [(f"C{i}", "C") for i in range(6)]),
            ("PHE", [("CG", "C"), ("CD1", "C"), ("CD2", "C"),
                     ("CE1", "C"), ("CE2", "C"), ("CZ", "C")]),
        ])
        traj = _trajectory(top, _ring(4.94) + _ring(0.12))
        chemistry = _chemistry("c1ccccc1")

        across = rules.pi_stacking(traj, chemistry, list(range(6)),
                                   list(range(6, 12)))
        ignored = rules.pi_stacking(traj, chemistry, list(range(6)),
                                    list(range(6, 12)), periodic=False)

        assert len(across) >= 1
        assert across[0].distance_nm == pytest.approx(0.18, abs=0.03)
        assert ignored == []


class TestAGroupSplitAcrossAFaceHasACentre:
    """Before any distance, the group itself has to be whole."""

    def test_a_carboxylate_straddling_the_face(self) -> None:
        """Two oxygens 0.22 nm apart, stored either side of x = 0. Their raw
        mean is the middle of the box, 2.5 nm from both of them."""
        top = _topology([("ASP", [("OD1", "O"), ("OD2", "O")])])
        traj = _trajectory(top, [[0.02, 2.5, 2.5], [4.89, 2.5, 2.5]])

        whole = rules._group_centres(traj, [[0, 1]], True)[0, 0]
        raw = rules._group_centres(traj, [[0, 1]], False)[0, 0]

        # The true centre is 0.09 nm from one oxygen, wrapping through zero.
        assert min(whole[0] % BOX_NM, BOX_NM - whole[0] % BOX_NM) < 0.12
        assert raw[0] == pytest.approx(2.455, abs=0.01)

    def test_a_ring_straddling_the_face_keeps_its_radius(self) -> None:
        top = _topology([("PHE", [(f"C{i}", "C") for i in range(6)])])
        points = _ring(0.0)  # centred on the face, so three atoms wrap
        wrapped = [[x % BOX_NM, y, z] for x, y, z in points]
        traj = _trajectory(top, wrapped)

        centres, normals = rules._ring_geometry(traj, [list(range(6))], True)
        rebuilt = rules._ring_geometry(traj, [list(range(6))], False)[0]

        # Every atom sits one ring radius from the whole centre.
        offsets = md.compute_displacements(
            traj, np.array([[0, i] for i in range(6)]), periodic=True)[0]
        spread = np.linalg.norm(
            traj.xyz[0, 0] + offsets - centres[0, 0], axis=-1)
        assert np.allclose(spread, 0.14, atol=0.01)
        assert np.isclose(np.linalg.norm(normals[0, 0]), 1.0)

        # And the raw centre is nowhere near the ring.
        assert abs(rebuilt[0, 0][0] - 2.5) < 1.0


class TestTheTrajectoryIsMadeWholeWhenItIsLoaded:
    """Nothing downstream images anything, so it has to happen at load."""

    def _split_chain(self):
        """Four bonded carbons 1 A apart, stored either side of x = 0."""
        top = md.Topology()
        residue = top.add_residue("LIG", top.add_chain())
        atoms = [top.add_atom(f"C{i}", md.element.carbon, residue)
                 for i in range(4)]
        for first, second in zip(atoms, atoms[1:]):
            top.add_bond(first, second)

        xyz = np.array([[[4.85, 2.5, 2.5], [4.95, 2.5, 2.5],
                         [0.05, 2.5, 2.5], [0.15, 2.5, 2.5]]],
                       dtype=np.float32)
        traj = md.Trajectory(xyz, top)
        traj.unitcell_vectors = np.array([np.eye(3) * BOX_NM],
                                         dtype=np.float32)
        return traj

    def test_a_split_molecule_is_made_whole(self) -> None:
        from fastmdxplora.analysis.loading import _made_whole

        whole = _made_whole(self._split_chain())

        extent = float(np.ptp(whole.xyz[0][:, 0]))
        assert extent == pytest.approx(0.30, abs=0.05), (
            f"four atoms 1 A apart span {extent:.2f} nm, so the molecule is "
            f"still split across the box face"
        )

    def test_the_radius_of_gyration_is_the_molecule_and_not_the_box(
        self,
    ) -> None:
        """The measure that read 2.643 nm for a 0.723 nm globule."""
        from fastmdxplora.analysis.loading import _made_whole

        split = self._split_chain()
        before = float(md.compute_rg(split)[0])
        after = float(md.compute_rg(_made_whole(split))[0])

        assert before > 2.0     # what every analysis saw
        assert after < 0.2      # the molecule it actually is

    def test_it_reaches_the_loader(self, tmp_path) -> None:
        """End to end, on a residue whose bonds mdtraj infers from the PDB."""
        from fastmdxplora.analysis.loading import load_trajectory

        top = md.Topology()
        residue = top.add_residue("ALA", top.add_chain())
        for name, element in (("N", "N"), ("CA", "C"), ("CB", "C"),
                              ("C", "C"), ("O", "O")):
            top.add_atom(name, md.element.get_by_symbol(element), residue)
        # Spread either side of the face, one bond length apart.
        xyz = np.array([[[4.85, 2.5, 2.5], [4.95, 2.5, 2.5], [0.05, 2.5, 2.5],
                         [0.15, 2.5, 2.5], [0.25, 2.5, 2.5]]],
                       dtype=np.float32)
        traj = md.Trajectory(xyz, top)
        traj.unitcell_vectors = np.array([np.eye(3) * BOX_NM],
                                         dtype=np.float32)
        path = tmp_path / "split.pdb"
        traj.save_pdb(str(path))

        loaded = load_trajectory(path)

        assert float(np.ptp(loaded.xyz[0][:, 0])) < 0.6

    def test_a_topology_without_bonds_says_so_rather_than_guessing(
        self, tmp_path, caplog
    ) -> None:
        """A molecule cannot be made whole without knowing what a molecule
        is. Loading still succeeds -- an analysis of what was loaded beats
        refusing to load it -- and the reason is on the record, so a run that
        went unimaged can be recognised afterwards."""
        import logging

        from fastmdxplora.analysis.loading import _made_whole

        top = md.Topology()
        residue = top.add_residue("LIG", top.add_chain())
        for i in range(3):
            top.add_atom(f"C{i}", md.element.carbon, residue)
        traj = md.Trajectory(
            np.array([[[4.9, 2.5, 2.5], [0.05, 2.5, 2.5], [0.1, 2.5, 2.5]]],
                     dtype=np.float32), top)
        traj.unitcell_vectors = np.array([np.eye(3) * BOX_NM],
                                         dtype=np.float32)

        # Attached to the logger that emits, not to root: `setup_console`
        # sets propagate = False on `fastmdx`, and the TODO records this
        # exact mistake being made three times.
        logger = logging.getLogger("fastmdx.analysis.loading")
        with caplog.at_level(logging.WARNING, logger=logger.name):
            logger.addHandler(caplog.handler)
            try:
                result = _made_whole(traj)
            finally:
                logger.removeHandler(caplog.handler)

        assert result is not None
        assert "periodic boundary" in caplog.text

    def test_a_trajectory_with_no_box_is_left_alone(self, tmp_path) -> None:
        top = md.Topology()
        residue = top.add_residue("ALA", top.add_chain())
        for i in range(3):
            top.add_atom(f"C{i}", md.element.carbon, residue)
        xyz = np.array([[[0.0, 0.0, 0.0], [0.1, 0.0, 0.0], [0.2, 0.0, 0.0]]],
                       dtype=np.float32)
        path = tmp_path / "vacuum.pdb"
        md.Trajectory(xyz, top).save_pdb(str(path))

        from fastmdxplora.analysis.loading import load_trajectory

        loaded = load_trajectory(path)
        assert np.allclose(loaded.xyz[0], xyz[0], atol=1e-5)


class TestAWaterSiteFoundThroughAFaceIsOneSite:
    """AUD8: accepted on a minimum-image distance, recorded raw.

    A water sitting 0.28 nm from a site atom through the box face was found
    every frame and then stored wherever the file happened to put it. When
    the file alternated images, DBSCAN received two clusters a box length
    apart: one physical site reported as two, each at half occupancy, one of
    them at a coordinate outside the box. Below the default
    `minimum_occupancy` the smaller half is dropped and the survivor keeps
    the halved number.
    """

    def _run(self, tmp_path, alternate: bool):
        top = md.Topology()
        chain = top.add_chain()
        protein = top.add_residue("ALA", chain)
        top.add_atom("CA", md.element.carbon, protein)
        water = top.add_residue("HOH", chain)
        top.add_atom("O", md.element.oxygen, water)

        frames = []
        for i in range(20):
            # The water is always ~0.15 nm from the site through the face.
            here = 0.05 if (not alternate or i % 2 == 0) else 0.05 + BOX_NM
            frames.append([[4.90, 2.5, 2.5], [here % BOX_NM, 2.5, 2.5]])
        traj = md.Trajectory(np.array(frames, dtype=np.float32), top)
        traj.unitcell_vectors = np.repeat(
            np.array([np.eye(3) * BOX_NM], dtype=np.float32), 20, axis=0)
        return traj

    def test_the_positions_recorded_are_all_beside_the_site(
        self, tmp_path
    ) -> None:
        """Whatever image the file stored them in."""
        from fastmdxplora.analysis.water_sites import WaterSites

        traj = self._run(tmp_path, alternate=True)
        analysis = WaterSites(cutoff_nm=0.35, minimum_occupancy=0.2)
        analysis.site_selection = "name CA"
        result = analysis.compute(traj)

        assert len(result) == 1, (
            f"one water at one site gave {len(result)} sites; an image swap "
            f"between frames split it"
        )
        row = result.iloc[0]
        assert float(row["occupancy"]) == pytest.approx(1.0, abs=0.05)

        # And the site sits beside the site atom rather than across the box.
        # Recorded contiguous with the protein, as an imaged trajectory is,
        # so the coordinate may lie just outside [0, L) -- what matters is
        # that it is 1.5 A from the atom it is a site for, not 48.5.
        site_atom = np.array([4.90, 2.5, 2.5])
        here = np.array([float(row["x"]), float(row["y"]), float(row["z"])])
        assert float(np.linalg.norm(here - site_atom)) < 0.35
