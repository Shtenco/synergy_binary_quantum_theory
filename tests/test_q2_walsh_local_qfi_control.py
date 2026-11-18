from scripts.q2_walsh_local_qfi_control import run


def test_q2_walsh_local_face_qfi_isotropic():
    result = run()
    assert result["passed"] is True
    import numpy as np
    assert np.allclose(
        np.asarray(result["four_face_QFI_sum"]),
        (8.0 / 3.0) * np.eye(3),
        atol=1e-12,
    )


def test_global_su2_gauge_qfi_vanishes_after_gauss_projection():
    result = run()
    import numpy as np
    assert result["gauss_singlet_weight"] == pytest.approx(2.0 / 9.0, abs=1e-12)
    assert result["global_SU2_generator_residual"] < 1e-12
    assert np.linalg.norm(result["global_SU2_QFI_after_Gauss_projection"]) < 1e-12


def test_relative_product_probe_is_only_a_uniform_local_control():
    result = run()
    assert len(result["relative_product_probe"]) == 6
    assert all(
        row["reduced_probe_conductance_FQ_over_4"] == pytest.approx(1.0 / 3.0, abs=1e-12)
        for row in result["relative_product_probe"]
    )
    assert "not a BQG QFI conductance" in result["claim_boundary"]


import pytest
