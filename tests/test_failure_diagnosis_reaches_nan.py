"""A NaN raised by the integrator reaches the diagnosis.

The machinery to read a failed state existed and was wired to the check that
runs at a stage boundary. OpenMM detects a non-finite coordinate during
integration and throws, so that check never ran -- and that is the common way
a simulation dies, not the rare one. A real run hit it and got the generic
list of settings to try, which is what the diagnosis was written to replace.
"""
from __future__ import annotations
import math, types, pytest
from fastmdxplora.simulation.runner import _state_for_diagnosis, _validation_error


class _Res:
    def __init__(self, name, idx): self.name, self.index = name, idx
class _Atom:
    def __init__(self, res, idx): self.residue, self.index = res, idx
class _Topo:
    def __init__(self, atoms): self._a = atoms
    def atoms(self): return iter(self._a)
    def getNumAtoms(self): return len(self._a)


def _topology(n=10, bad=3):
    res = _Res("TRP", 0)
    return _Topo([_Atom(res, i) for i in range(n)])


def _positions(n=10, bad=3):
    return [[float("nan")] * 3 if i < bad else [0.1 * i] * 3 for i in range(n)]


class _Ctx:
    def __init__(self, ok=True): self.ok = ok
    def getState(self, **kw):
        if not self.ok: raise RuntimeError("context gone")
        return types.SimpleNamespace(getPositions=lambda **k: _positions())


class _Sim:
    def __init__(self, ok=True):
        self.context, self.topology = _Ctx(ok), _topology()


class TestTheFailedStateIsRecovered:
    def test_positions_come_back_from_a_live_context(self):
        topo, pos = _state_for_diagnosis(_Sim())
        assert topo is not None and pos is not None
        assert any(math.isnan(v) for row in pos for v in row)

    def test_an_unreadable_context_is_not_an_error(self):
        """A diagnosis that fails must not replace the failure being reported."""
        assert _state_for_diagnosis(_Sim(ok=False)) == (None, None)


class TestTheMessageCarriesBoth:
    """`diagnose_failure` needs OpenMM to read a topology, and where it is
    absent the general message is the right fallback. What is tested here is
    the wiring -- that a state reaches the diagnosis at all -- because that
    is what was missing, and it must hold whether or not OpenMM is
    importable in the environment running the tests."""

    def test_the_state_is_handed_over(self, monkeypatch) -> None:
        seen = {}

        # `**rest` because the real one grew a `platform` argument -- the
        # advice it offers depends on which platform is running, since there
        # is no `Precision` property to change on CPU. A double that refuses
        # an argument the real function takes raises TypeError inside a
        # `try`, and the caller falls through to its generic message: the
        # diagnosis silently stops being reached, which is exactly the fault
        # this file exists to catch.
        def _fake(topology, positions, *, stage, **rest):
            seen["topology"], seen["stage"] = topology, stage
            return types.SimpleNamespace(
                as_text=lambda: f"The simulation became unstable during {stage}. TRP")

        import fastmdxplora.simulation.diagnose as d
        monkeypatch.setattr(d, "diagnose_failure", _fake)
        err = _validation_error("NPT equilibration", "Particle coordinate is NaN",
                                topology=_topology(), positions=_positions())
        assert seen["stage"] == "NPT equilibration"
        assert "unstable during NPT equilibration" in str(err)

    def test_openmm_own_words_survive(self, monkeypatch) -> None:
        """Searchable, and they link to its FAQ. The diagnosis says which
        atoms; neither substitutes for the other."""
        import fastmdxplora.simulation.diagnose as d
        monkeypatch.setattr(
            d, "diagnose_failure",
            lambda t, p, *, stage, **rest: types.SimpleNamespace(
                as_text=lambda: "DIAG"))
        err = _validation_error("NPT equilibration", "Particle coordinate is NaN",
                                topology=_topology(), positions=_positions())
        assert "DIAG" in str(err)
        assert "OpenMM reported: Particle coordinate is NaN" in str(err)

    def test_a_diagnosis_that_fails_does_not_replace_the_failure(
            self, monkeypatch) -> None:
        import fastmdxplora.simulation.diagnose as d

        def _boom(*a, **k): raise RuntimeError("no openmm")
        monkeypatch.setattr(d, "diagnose_failure", _boom)
        err = _validation_error("NPT equilibration", "Particle coordinate is NaN",
                                topology=_topology(), positions=_positions())
        assert "Try safer settings" in str(err)

    def test_without_a_state_the_general_message_still_appears(self) -> None:
        err = _validation_error("NPT equilibration", "boom")
        assert "Try safer settings" in str(err)
