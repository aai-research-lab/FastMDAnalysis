"""Runs that differ only in how they are driven still solvate separately.

`explorer` shares one prepared system across umbrella windows, for a reason
its own docstring states: solvation does not place water the same way
twice, and water arranged differently between windows is noise rather than
physics. That path is gated on `_is_umbrella`.

A seed sweep is the same case in kind -- one system, one setup block, only
the integrator seed differing -- and is not umbrella, so each member
solvates independently. Measured on a three-member sweep: 30,654 atoms
against 30,803. Anything called a replica sweep then measures dynamics
variance *plus* solvation variance, and no recorded setting shows it,
because the difference is in water that no setting names.

The software says so at plan time. It does not share setup by default:
that would silently change what an existing config produces, and a study
half-run under the old behaviour would stop being comparable with its own
earlier members.
"""

from __future__ import annotations

import logging

import pytest

from fastmdxplora.batch.explorer import (
    _say_if_the_replicas_will_not_share_water,
)


class _Spec:
    """The two fields the notice reads, and nothing else."""

    def __init__(self, system="1L2Y.pdb", setup=None, simulation=None):
        self.system = system
        self.options = {"setup": setup or {}, "simulation": simulation or {}}


def _warnings(caplog, specs):
    caplog.clear()
    with caplog.at_level(logging.WARNING, logger="fastmdxplora.batch.explorer"):
        _say_if_the_replicas_will_not_share_water(specs)
    return [r.getMessage() for r in caplog.records]


class TestASeedSweepIsTold:
    def test_it_names_the_setting(self, caplog):
        specs = [_Spec(simulation={"random_seed": s}) for s in (1, 2, 3)]
        messages = _warnings(caplog, specs)

        assert len(messages) == 1
        assert "prepared_from" in messages[0]
        assert "3 runs" in messages[0]

    def test_silence_once_prepared_from_is_set(self, caplog):
        """The remedy having been applied, there is nothing to say."""
        specs = [_Spec(simulation={"random_seed": s, "prepared_from": "out/setup"})
                 for s in (1, 2, 3)]

        assert _warnings(caplog, specs) == []

    def test_different_systems_are_different_systems(self, caplog):
        """Not a replica set. Sharing one preparation would be wrong here,
        so saying they do not share it would be noise."""
        specs = [_Spec(system=s) for s in ("1L2Y.pdb", "1UBQ.pdb", "3PTB.pdb")]

        assert _warnings(caplog, specs) == []

    def test_a_sweep_over_setup_is_left_alone(self, caplog):
        """Swept preparation settings mean the runs are *meant* to be
        prepared differently."""
        specs = [_Spec(setup={"padding_nm": p}) for p in (1.0, 1.2, 1.4)]

        assert _warnings(caplog, specs) == []

    def test_one_run_is_not_a_replica_set(self, caplog):
        assert _warnings(caplog, [_Spec()]) == []


class TestItDoesNotChangeWhatRuns:
    def test_nothing_is_written_into_the_specs(self, caplog):
        """Saying so must not become doing so: a config's meaning is
        unchanged by being warned about."""
        specs = [_Spec(simulation={"random_seed": s}) for s in (1, 2)]
        before = [dict(s.options["simulation"]) for s in specs]

        _warnings(caplog, specs)

        assert [dict(s.options["simulation"]) for s in specs] == before
