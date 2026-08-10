"""Advice that names a setting the platform does not have costs a run.

`Precision` is a property of OpenMM's GPU platforms. The CPU platform offers
`Threads` and `DeterministicForces` and nothing else, so `precision=mixed`
there is a request nothing acts on -- and the banner said "CPU
(precision=mixed)" as though it had.

That matters most where it is least affordable. After a non-finite coordinate,
the diagnosis offers remedies to try, and one of them was "run in double
precision". On CPU that changes nothing: somebody chasing an unstable run
would spend another hour and learn nothing about their system.
"""

from __future__ import annotations

import numpy as np
import pytest

from fastmdxplora.simulation.diagnose import diagnose_failure


openmm_app = pytest.importorskip("openmm.app")


def _water_box(n: int = 60):
    """A real OpenMM topology: the diagnosis type-checks for one."""
    topology = openmm_app.Topology()
    chain = topology.addChain()
    for _ in range(n):
        residue = topology.addResidue("HOH", chain)
        topology.addAtom("O", openmm_app.element.oxygen, residue)
    return topology


def _positions(n: int = 60, bad: int = 55) -> np.ndarray:
    """Most of the box non-finite: the branch that offers general advice."""
    xyz = np.ones((n, 3), dtype=float)
    xyz[:bad] = np.nan
    return xyz


class TestThePrecisionAdviceIsOfferedWhereItCanBeTaken:
    def test_not_on_cpu(self) -> None:
        text = diagnose_failure(_water_box(), _positions(), stage="Production",
                                platform="CPU").as_text()
        assert "--simulate-precision" not in text

    @pytest.mark.parametrize("platform", ["CUDA", "OpenCL", "HIP"])
    def test_yes_on_a_gpu(self, platform: str) -> None:
        text = diagnose_failure(_water_box(), _positions(), stage="Production",
                                platform=platform).as_text()
        assert "--simulate-precision double" in text

    def test_an_unknown_platform_keeps_it(self) -> None:
        """Withholding advice because the platform could not be identified
        would lose it on every platform this does not know about."""
        text = diagnose_failure(_water_box(), _positions(),
                                stage="Production").as_text()
        assert "--simulate-precision double" in text

    def test_the_other_advice_survives(self) -> None:
        """Only the one remedy that cannot be taken is withheld."""
        text = diagnose_failure(_water_box(), _positions(), stage="Production",
                                platform="CPU").as_text()
        assert "--simulate-timestep-fs" in text
        assert "--simulate-friction-per-ps" in text


class TestTheBannerSaysWhatWasApplied:
    def test_it_reports_the_applied_precision_not_the_request(self) -> None:
        from pathlib import Path
        import fastmdxplora.simulation.runner as runner

        source = Path(runner.__file__).read_text(encoding="utf-8")
        assert 'applied = props.get("Precision")' in source
        assert "this platform has no such " in source

    def test_the_platform_is_read_from_the_context(self) -> None:
        """Rather than threaded down from where it was chosen: the context is
        the thing that knows what it is running on."""
        from pathlib import Path
        import fastmdxplora.simulation.runner as runner

        source = Path(runner.__file__).read_text(encoding="utf-8")
        assert "simulation.context.getPlatform().getName()" in source
