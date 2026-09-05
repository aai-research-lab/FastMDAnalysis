"""A platform that loads is not a platform that runs.

`auto` tries CUDA, then OpenCL, then CPU, and it probes each GPU platform
before choosing it. The probe built a Context around a System holding one
particle and no forces. That compiles almost nothing, so a platform whose
kernels will not build passed the probe and failed on the real system a
second later.

Seen on a workstation whose driver rejected the installed OpenMM build's
PTX. `openmm.testInstallation` reported `CUDA - Error computing forces`,
while `auto` selected CUDA, and every run died with
`CUDA_ERROR_UNSUPPORTED_PTX_VERSION` -- with a working OpenCL platform
sitting next in the candidate list, never tried.

The probe now computes a force, which is what `testInstallation` does and
what distinguishes the two outcomes. These tests pin that: a platform that
raises only once forces are requested must be rejected, and the next
candidate must be chosen instead.
"""

from __future__ import annotations


from fastmdxplora.simulation import runner


class _FakeState:
    def __init__(self, fail):
        self._fail = fail

    def getForces(self, asNumpy=False):
        if self._fail:
            raise RuntimeError(
                "Error loading CUDA module: CUDA_ERROR_UNSUPPORTED_PTX_VERSION (222)")
        return [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]


class _FakeContext:
    """Builds fine. Fails when asked for forces, like a bad PTX build."""

    def __init__(self, fail_on_forces):
        self._fail = fail_on_forces

    def setPositions(self, positions):
        return None

    def getState(self, getForces=False):
        return _FakeState(self._fail)


def _fake_openmm(broken: set[str]):
    """An OpenMM whose named platforms fail only when forces are computed."""

    class _Force:
        def addParticle(self, *a):
            return None

    class _System:
        def addParticle(self, *a):
            return None

        def addForce(self, *a):
            return None

    def _context(system, integrator, platform, props):
        return _FakeContext(platform.getName() in broken)

    class _MM:
        System = _System
        NonbondedForce = _Force
        Context = staticmethod(_context)

        @staticmethod
        def VerletIntegrator(step):
            return object()

    class _Quantity:
        """Real OpenMM units support `list * unit`; a bare float does not.

        The first version of this fake used floats, and `_probe_platform`
        rejected a perfectly good platform with
        "can't multiply sequence by non-int of type 'float'" -- a fault in
        the fixture that would have read as a fault in the code.
        """

        def __rmul__(self, other):
            return other

        def __mul__(self, other):
            return other

    class _Unit:
        picoseconds = _Quantity()
        nanometers = _Quantity()

    return {"openmm": _MM(), "unit": _Unit()}


class TestTheProbeComputesAForce:
    def test_a_platform_that_fails_on_forces_is_rejected(self, monkeypatch):
        """The exact shape of the workstation failure: the Context builds,
        the kernels do not compile until forces are asked for."""
        omm = _fake_openmm(broken={"CUDA"})
        platform = type("P", (), {"getName": lambda self: "CUDA"})()

        monkeypatch.setattr(runner, "omm", omm, raising=False)
        usable = runner._probe_platform(omm, platform, "CUDA", {})

        assert usable is False

    def test_a_working_platform_is_accepted(self, monkeypatch):
        omm = _fake_openmm(broken=set())
        platform = type("P", (), {"getName": lambda self: "OpenCL"})()

        assert runner._probe_platform(omm, platform, "OpenCL", {}) is True

    def test_cpu_is_not_probed(self):
        """CPU and Reference are always usable and cost nothing to trust."""
        assert runner._probe_platform(None, None, "CPU", {}) is True
        assert runner._probe_platform(None, None, "Reference", {}) is True
