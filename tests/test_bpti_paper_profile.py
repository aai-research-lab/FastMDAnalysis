"""Regression tests for the complete BPTI paper reproduction profile."""

from __future__ import annotations

from pathlib import Path

import yaml

from fastmdxplora.config import apply_profile, validate_config


def test_bpti_paper_profile_has_the_expected_phase_settings():
    profile = apply_profile({"preset": "bpti-paper"})
    validate_config(profile, require_systems=True)

    assert profile["systems"] == [{"id": "bpti-6pti", "system": "6PTI"}]
    assert profile["setup"]["forcefield"] == "charmm36m"
    assert profile["setup"]["hydrogen_mass_amu"] == 1.5
    assert profile["simulation"]["nvt_steps"] == 250000
    assert profile["simulation"]["npt_steps"] == 1000000
    assert profile["simulation"]["trajectory_interval_steps"] == 5000
    assert profile["analysis"]["stride"] == 2
    assert profile["analysis"]["options"]["cluster"] == {
        "methods": ["hierarchical"], "n_clusters": 6
    }
    assert profile["analysis"]["options"]["dimred"] == {"methods": ["pca"]}
    assert profile["analysis"]["options"]["rmsd"]["selection"] == "protein and backbone"
    assert profile["analysis"]["options"]["rmsf"]["align"] is False


def test_bpti_example_is_valid_and_self_contained():
    root = Path(__file__).resolve().parents[1]
    data = yaml.safe_load((root / "examples" / "bpti_paper_6pti.yaml").read_text())
    expanded = apply_profile(data)
    validate_config(expanded, require_systems=True)
    assert expanded["simulation"]["production_steps"] == 50000000
