import json

import pytest

from scripts.bqg_microscopic_node_refinement_gate import run as run_fixed_pure_ancilla
from scripts.bqg_microscopic_bilinear_refinement_gate import run as run_bilinear_refinement


def test_fixed_pure_ancilla_shortcut_remains_an_explicit_scientific_no_go():
    result = run_fixed_pure_ancilla()

    # The no-go is the expected scientific result, not a test-suite failure.
    assert result["passed"] is False
    assert result["relative_reconstruction_error"] == pytest.approx(1.0, abs=1e-10)
    assert result["phase_fixed_error"] == pytest.approx(2**0.5, abs=1e-10)
    assert result["Gauss_projection_leakage_sum"] < 1e-9
    assert "fixed pure" in result["claim_boundary"].lower()


def test_bilinear_ancilla_refinement_is_an_s4_intertwiner(capsys):
    exit_code = run_bilinear_refinement()
    result = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert result["passed"] is True
    assert result["coarse_22_tensor_rank"] == 2
    assert result["ancilla_invariant_vector_dimension"] == 0
    assert result["S4_covariance_error"] < 1e-9
    assert result["fusion_source_projector_error"] < 1e-9
    assert result["fusion_target_projector_error"] < 1e-9
    assert result["Gauss_projection_leakage"] < 1e-9
