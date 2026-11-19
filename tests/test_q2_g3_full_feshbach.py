import numpy as np

from scripts.q2_g3_full_feshbach import full_feshbach_matrices, matrix_diagnostics
from scripts.q2_matrix_free_feshbach import MatrixFreeFeshbach


def _dense_reference(model):
    n, m = model.n_vertices, model.n_fibres
    u = np.column_stack([model.apply_u(np.eye(m)[:, j]) for j in range(m)])
    laplacian = model.laplacian.toarray()
    p = u @ u.T
    q = np.eye(n) - p
    b = q @ laplacian @ u
    h = q @ laplacian @ q + p
    inverse = np.linalg.inv(h)
    return b.T @ inverse @ b, b.T @ inverse @ inverse @ b


def test_full_feshbach_matrices_match_dense_small_graph():
    edges = np.array([[0, 1], [1, 2], [2, 3], [3, 0]], dtype=np.uint32)
    labels = np.array([0, 0, 1, 1], dtype=np.uint16)
    model = MatrixFreeFeshbach(edges, labels)
    expected_correction, expected_slope = _dense_reference(model)

    correction, slope, info = full_feshbach_matrices(
        model, rtol=1e-12, atol=1e-13, maxiter=100
    )

    assert np.allclose(correction, expected_correction, atol=1e-10)
    assert np.allclose(slope, expected_slope, atol=1e-10)
    assert info["columns_solved"] == model.n_fibres
    assert info["max_first_solve_residual"] < 1e-10
    assert info["max_second_solve_residual"] < 1e-10
    assert np.allclose(correction, correction.T, atol=1e-10)
    assert np.allclose(slope, slope.T, atol=1e-10)


def test_full_feshbach_matrix_diagnostics_reject_nonfinite_and_nonsquare():
    import pytest

    with pytest.raises(ValueError, match="square"):
        matrix_diagnostics(np.ones((2, 3)))
    with pytest.raises(ValueError, match="non-finite"):
        matrix_diagnostics(np.array([[np.nan]]))


def test_full_feshbach_matrix_diagnostics_report_psd_spectrum():
    report = matrix_diagnostics(np.array([[1.0, -1.0], [-1.0, 1.0]]))
    assert report["symmetry_error"] == 0.0
    assert report["minimum_eigenvalue"] >= -1e-14
    assert report["maximum_eigenvalue"] == 2.0
    assert report["trace"] == 2.0
