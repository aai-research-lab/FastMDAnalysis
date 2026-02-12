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
    assert analysis.data.shape[0] == 1
    assert "phi_residues" in analysis.results
    assert list(np.asarray(analysis.results["phi_residues"]).astype(int)) == [0]
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
    # Ensure the combined analysis truly computed only 2 residues worth of data
    assert analysis.results["phi_avg"].shape[0] == 2
    assert analysis.results["psi_avg"].shape[0] == 2
    assert analysis.results["omega_avg"].shape[0] == 2
    assert list(np.asarray(analysis.results["phi_residues"]).astype(int)) == [0, 1]
    assert list(np.asarray(analysis.results["psi_residues"]).astype(int)) == [0, 1]
    assert list(np.asarray(analysis.results["omega_residues"]).astype(int)) == [0, 1]
    for key in ("phi_avg_filtered", "psi_avg_filtered", "omega_avg_filtered"):
        assert key in analysis.results
        filtered = analysis.results[key]
        assert filtered.shape[0] == 2
        base_key = key.replace("_filtered", "")
        np.testing.assert_allclose(filtered[:, 0], analysis.results[base_key][:2, 0], atol=1e-6)
