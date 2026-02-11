import warnings
import numpy as np


def _unknown_warnings(messages):
    return [m for m in messages if "Unknown options" in m or "Unsupported options" in m]


def test_phi_residue_alias_no_unknown_warning(fastmda):
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        analysis = fastmda.phi(residue=0)

    messages = [str(wi.message) for wi in w]
    assert not _unknown_warnings(messages)

    assert "phi_avg_filtered" in analysis.results
    filtered = analysis.results["phi_avg_filtered"]
    assert filtered.shape[0] == 1
    np.testing.assert_allclose(filtered[0], analysis.data[0], atol=1e-6)


def test_dihedrals_residue_selection_alias_propagates(fastmda):
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        analysis = fastmda.dihedrals(residue_selection=[0, 1])

    messages = [str(wi.message) for wi in w]
    assert not _unknown_warnings(messages)

    for key in ("phi_avg_filtered", "psi_avg_filtered", "omega_avg_filtered"):
        assert key in analysis.results
        filtered = analysis.results[key]
        assert filtered.shape[0] == 2
        base_key = key.replace("_filtered", "")
        np.testing.assert_allclose(filtered[:, 0], analysis.results[base_key][:2, 0], atol=1e-6)
