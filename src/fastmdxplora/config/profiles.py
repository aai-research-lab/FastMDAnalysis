"""Named whole-workflow configuration profiles."""

from __future__ import annotations

from copy import deepcopy
from typing import Any


PROFILES: dict[str, dict[str, Any]] = {
    "bpti-paper": {
        "systems": [{"id": "bpti-6pti", "system": "6PTI"}],
        "setup": {
            "ph": 7.0,
            "keep_heterogens": False,
            "keep_water": False,
            "forcefield": "charmm36m",
            "water_model": "tip3p",
            "solvent_padding_nm": 1.0,
            "box_shape": "cube",
            "ion_positive": "Na+",
            "ion_negative": "Cl-",
            "ion_concentration_M": 0.15,
            "nonbonded_method": "PME",
            # OpenMM exposes one shared cutoff. The paper's asymmetric
            # 1.0 nm electrostatic / 1.2 nm LJ scheme is documented in the
            # companion guide; do not pretend a zero-width switch is valid.
            "nonbonded_cutoff_nm": 1.0,
            "use_switching_function": False,
            "switch_distance_nm": None,
            "constraints": "HBonds",
            "hydrogen_mass_amu": 1.5,
            "temperature_K": 300.0,
        },
        "simulation": {
            "duration_ns": 100.0,
            "timestep_fs": 2.0,
            "temperature_K": 300.0,
            "friction_per_ps": 1.0,
            "nvt_steps": 250000,
            "npt_steps": 1000000,
            "production_steps": 50000000,
            "barostat_frequency": 50,
            "platform": "CUDA",
            "precision": "mixed",
            "trajectory_interval_steps": 5000,
            "checkpoint_interval_steps": 100000,
            "live_telemetry": False,
        },
        "analysis": {
            "scope": "protein",
            "stride": 2,
            "include": ["rmsd", "rmsf", "rg", "hbonds", "sasa", "ss", "cluster", "dimred"],
            "options": {
                "rmsd": {"selection": "protein and backbone", "ref": 0, "align": True},
                "rmsf": {"selection": "protein", "per_residue": False, "align": False},
                "sasa": {"mode": "total", "probe_radius": 0.14},
                "cluster": {"methods": ["hierarchical"], "n_clusters": 6},
                "dimred": {"methods": ["pca"]},
            },
        },
    },
}


def available_profiles() -> tuple[str, ...]:
    """Return profile names in stable order."""
    return tuple(sorted(PROFILES))


def _merge(base: dict[str, Any], override: dict[str, Any]) -> dict[str, Any]:
    """Recursively merge mappings; scalar values and lists are replaced."""
    merged = deepcopy(base)
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(merged.get(key), dict):
            merged[key] = _merge(merged[key], value)
        else:
            merged[key] = deepcopy(value)
    return merged


def apply_profile(data: dict[str, Any]) -> dict[str, Any]:
    """Expand ``preset:`` while allowing explicit YAML/CLI values to win."""
    preset = data.get("preset")
    if preset is None:
        return deepcopy(data)
    key = str(preset).lower()
    if key not in PROFILES:
        valid = ", ".join(available_profiles())
        raise ValueError(f"Unknown project preset {preset!r}. Valid presets: {valid}.")
    explicit = deepcopy(data)
    explicit["preset"] = key
    return _merge(PROFILES[key], explicit)
