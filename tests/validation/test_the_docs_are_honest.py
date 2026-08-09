"""The documentation says things about the code that the code can check.

Four claims in the published docs had gone stale without anything noticing:
the number of analyses, the number of collective variables in three separate
places, which analyses need a ligand, and which section three enhanced
sampling analyses belonged under. Each was correct when written. Each became
wrong when something was added, and a reader had no way to tell.

Prose has to be written by hand. A count does not, and a count is what goes
wrong, so the counts are pinned here.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

import fastmdxplora.analysis.analyze  # noqa: F401  -- registers the analyses
from fastmdxplora.analysis.orchestrator import _REGISTRY
from fastmdxplora.simulation.metadynamics import COLLECTIVE_VARIABLES

DOCS = Path(__file__).resolve().parents[2] / "docs"
ROOT = Path(__file__).resolve().parents[2]

#: The three that read a biased run's own result rather than the trajectory.
METHOD_GATED = ("pmf", "metad_surface", "steered_work")


def _analysis_section() -> str:
    text = (DOCS / "phases.md").read_text(encoding="utf-8")
    start = text.index("## analysis")
    return text[start:text.index("## report", start)]


class TestEveryAnalysisIsDocumented:
    def test_the_table_lists_exactly_what_is_registered(self) -> None:
        documented = set(re.findall(r"^\| `([a-z_]+)` \|", _analysis_section(),
                                    re.M))
        assert documented == set(_REGISTRY), (
            "docs/phases.md and the analysis registry disagree; "
            f"undocumented={sorted(set(_REGISTRY) - documented)}, "
            f"documented but absent={sorted(documented - set(_REGISTRY))}")

    def test_the_stated_count_is_the_real_one(self) -> None:
        """'Fifteen analyses' outlived three additions."""
        section = _analysis_section()
        trajectory = len(_REGISTRY) - len(METHOD_GATED)
        words = {9: "Nine", 10: "Ten", 14: "Fourteen", 15: "Fifteen",
                 16: "Sixteen", 17: "Seventeen", 18: "Eighteen",
                 19: "Nineteen", 20: "Twenty"}
        assert f"{words[trajectory]} analyses" in section, (
            f"there are {trajectory} analyses of the trajectory; "
            "the sentence introducing them says otherwise")

    def test_the_ligand_gated_count_is_right(self) -> None:
        needs_ligand = [n for n, c in _REGISTRY.items()
                        if getattr(c, "requires_ligand", False)]
        assert len(needs_ligand) == 5, (
            "the docs say five ligand analyses run automatically; "
            f"the registry has {len(needs_ligand)}: {sorted(needs_ligand)}")
        assert "five ligand analyses" in _analysis_section()


class TestTheEnhancedSamplingThreeAreNotLigandAnalyses:
    """They were filed under 'Protein and ligand together', which is where a
    row lands if it is appended to the end of the file. None of them needs a
    ligand and none runs on an ordinary trajectory."""

    def test_none_of_them_requires_a_ligand(self) -> None:
        for name in METHOD_GATED:
            assert not getattr(_REGISTRY[name], "requires_ligand", False), name

    def test_each_runs_only_after_its_own_method(self) -> None:
        gates = ("requires_umbrella", "requires_metadynamics",
                 "requires_steered")
        for name in METHOD_GATED:
            assert any(getattr(_REGISTRY[name], g, False) for g in gates), name

    def test_they_sit_under_their_own_heading(self) -> None:
        section = _analysis_section()
        heading = section.index("**Enhanced sampling**")
        ligand = section.index("**Protein and ligand together**")
        for name in METHOD_GATED:
            assert section.index(f"| `{name}` |") > heading, (
                f"{name} is documented above the Enhanced sampling heading")
        assert ligand < heading


class TestTheCollectiveVariablesAreAllNamed:
    """The count was wrong in three places at once, and the README disagreed
    with the page it links to."""

    @pytest.mark.parametrize("document", ["simulations.md", "phases.md"])
    def test_no_document_understates_them(self, document: str) -> None:
        text = (DOCS / document).read_text(encoding="utf-8")
        assert "five variables" not in text
        assert "on five variables" not in text

    def test_every_variable_is_named_where_they_are_listed(self) -> None:
        text = (DOCS / "simulations.md").read_text(encoding="utf-8")
        missing = [c for c in COLLECTIVE_VARIABLES if f"`{c}`" not in text]
        assert not missing, f"collective variables never named in docs: {missing}"

    def test_the_readme_and_the_page_agree(self) -> None:
        readme = (ROOT / "README.md").read_text(encoding="utf-8")
        words = {5: "five", 6: "six", 7: "seven", 8: "eight", 9: "nine",
                 10: "ten"}
        word = words[len(COLLECTIVE_VARIABLES)]
        assert word in readme.lower()
        assert word in (DOCS / "simulations.md").read_text(
            encoding="utf-8").lower()


class TestTheCountsInTheReadmeHold:
    def test_the_bilayers(self) -> None:
        from fastmdxplora.setup.membrane import LIPIDS
        assert len(LIPIDS) == 7, "the README says seven bilayers"
        text = (DOCS / "simulations.md").read_text(encoding="utf-8")
        for lipid in LIPIDS:
            assert lipid in text, f"{lipid} is offered and never documented"


class TestTheReweightingIsDescribedAsBuilt:
    """The estimator and its limits are claims about the code, and the code
    can check the ones that are structural."""

    @staticmethod
    def _section() -> str:
        """The section with its line wrapping flattened.

        Prose wraps where the column runs out, so "radius of gyration" spans
        a newline. A test that fails when a paragraph is rewrapped is a test
        that gets deleted rather than fixed."""
        import re
        text = (DOCS / "phases.md").read_text(encoding="utf-8")
        start = text.index("### Averages on a biased run")
        return re.sub(r"\s+", " ", text[start:text.index("## report", start)])

    def test_the_estimator_is_named_with_its_offset(self) -> None:
        """Without c(t) the weights rank frames by when they were written.
        Documenting exp(V/RT) alone would describe a different estimator."""
        section = self._section()
        assert "c(t)" in section
        assert "Tiwary" in section

    def test_every_reweighted_analysis_is_named(self) -> None:
        import re
        section = self._section().lower()
        for name, cls in _REGISTRY.items():
            if not getattr(cls, "reweightable", None):
                continue
            label = cls.reweightable[1].split(" (")[0].lower()
            assert label in section, (
                f"{name} is reweighted and the docs never say so")

    def test_what_is_not_corrected_is_said(self) -> None:
        section = self._section()
        assert "dimensionality reduction is not corrected" in section
        # Populations are corrected; which clusters exist is not.
        assert "which clusters exist" in section

    def test_the_two_methods_without_weights_are_explained(self) -> None:
        section = self._section()
        assert "umbrella window" in section.lower()
        assert "steered pull" in section.lower()
        assert "potential of mean force" in section

    def test_the_output_directory_is_documented(self) -> None:
        text = (DOCS / "phases.md").read_text(encoding="utf-8")
        assert "analysis/reweighted/" in text
        assert "reweighted_averages.json" in text
