"""How much of each window is discarded, and who decides it.

A window begins away from where it will settle, so the approach has to be
thrown out before its histogram is built. That was already happening -- a
fifth of every window -- but the fifth was hardcoded at the call site, absent
from the study record, and not a setting.

Which made it the one judgement in the umbrella study nobody could see or
change. `minimum_overlap` next to it is a setting for exactly this reason:
"a judgement about how much evidence a joint needs, so it belongs to whoever
is making the claim." How long a window takes to settle depends on the
barrier and the force constant, and is the same kind of judgement.
"""

from __future__ import annotations

import pytest

from fastmdxplora.simulation.umbrella import (
    _checked_fraction, collect_samples, plan_windows)


def _spec(**over):
    spec = {"from": 0.0, "to": 1.0, "n_windows": 3, "force_constant": 500.0,
            "collective_variable": "distance"}
    spec.update(over)
    return spec


class TestTheDiscardIsASetting:
    def test_a_fifth_by_default(self) -> None:
        assert plan_windows(_spec()).equilibration_fraction == pytest.approx(0.2)

    def test_a_study_can_choose_another(self) -> None:
        plan = plan_windows(_spec(equilibration_fraction=0.35))
        assert plan.equilibration_fraction == pytest.approx(0.35)

    def test_the_study_records_what_it_discarded(self) -> None:
        """A methods section states it, and a reader cannot recompute the
        free energy without knowing which samples went into it."""
        record = plan_windows(_spec(equilibration_fraction=0.3)).as_record()
        assert record["equilibration_fraction"] == pytest.approx(0.3)


class TestADegenerateFractionIsRefusedEarly:
    """Refused where the plan is made, not where the free energy comes out
    empty -- by then the windows have already run."""

    @pytest.mark.parametrize("bad", [0.0, 1.0, 1.5, -0.2])
    def test_it_must_leave_a_histogram_and_discard_something(self, bad) -> None:
        with pytest.raises(ValueError, match="equilibration_fraction"):
            _checked_fraction(bad)

    def test_zero_is_refused_rather_than_treated_as_none(self) -> None:
        """Keeping the approach to the centre pulls the histogram towards
        where the window started, which is the one place the free energy is
        guaranteed not to be flat."""
        with pytest.raises(ValueError, match="discarded"):
            _checked_fraction(0.0)

    def test_a_non_number_says_so(self) -> None:
        with pytest.raises(ValueError, match="must be a number"):
            _checked_fraction("a third")

    def test_the_plan_refuses_it_too(self) -> None:
        with pytest.raises(ValueError):
            plan_windows(_spec(equilibration_fraction=0.0))


class TestTheFractionReachesTheHistograms:
    def test_more_discarded_leaves_fewer_samples(self, tmp_path) -> None:
        window = tmp_path / "w0" / "simulation"
        window.mkdir(parents=True)
        (window / "COLVAR").write_text(
            "#! FIELDS time cv\n"
            + "\n".join(f"{i}.0 {i * 0.01:.4f}" for i in range(100)) + "\n",
            encoding="utf-8")

        kept = {f: collect_samples({0: tmp_path / "w0"},
                                   equilibration_fraction=f)[0].size
                for f in (0.1, 0.2, 0.5)}
        assert kept[0.1] == 90 and kept[0.2] == 80 and kept[0.5] == 50

    def test_the_call_site_passes_the_plans_fraction(self) -> None:
        """The default was hardcoded at the call site, so a study asking for
        a third would have reported a third and used a fifth."""
        from pathlib import Path
        import fastmdxplora.batch.explorer as explorer

        source = Path(explorer.__file__).read_text(encoding="utf-8")
        assert "equilibration_fraction=plan.equilibration_fraction" in source
