from scripts.q2_pl_generation3_preflight import barycentric_f_vector_3d, run


def test_generation3_exact_f_vector_and_dual_counts():
    result = run()
    assert result["passed"] is True
    assert result["generation3_f_vector"] == [40256, 261440, 442368, 221184]
    assert result["generation3_pl_dual_edges"] == 442368
    assert result["generation3_causal_active_edges"] == 512


def test_barycentric_f_vector_preserves_euler_characteristic_for_closed_3_complex():
    before = (1696, 10912, 18432, 9216)
    after = barycentric_f_vector_3d(before)
    assert before[0] - before[1] + before[2] - before[3] == 0
    assert after[0] - after[1] + after[2] - after[3] == 0


def test_preflight_marks_dense_spectral_cost_as_nontrivial():
    result = run()
    assert result["dense_spectral_matrix_raw_float64_bytes"] == 8 * 512 * 512
    assert result["jacobi_rotation_work_estimate_at_100_sweeps_fallback_only"] == 6697779200
    assert result["dense_lapack_cubic_work_proxy"] == 512 ** 3
    assert "numpy.linalg.eigvalsh" in result["large_spectrum_backend"]
    assert "No generation-3 complex" in result["claim_boundary"]

def test_barycentric_f_vector_rejects_wrong_dimension():
    import pytest

    with pytest.raises(ValueError, match="four f-vector entries"):
        barycentric_f_vector_3d((1, 2, 3))
