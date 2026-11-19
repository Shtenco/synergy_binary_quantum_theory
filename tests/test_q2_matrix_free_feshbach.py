import numpy as np
import pytest

scipy = pytest.importorskip("scipy")

from scripts.q2_matrix_free_feshbach import MatrixFreeFeshbach


@pytest.fixture
def cycle_operator():
    edges = np.array([[0, 1], [1, 2], [2, 3], [3, 0]], dtype=np.uint32)
    labels = np.array([0, 0, 1, 1], dtype=np.uint16)
    return MatrixFreeFeshbach(edges, labels)


def _dense_reference(operator):
    n = operator.n_vertices
    m = operator.n_fibres
    u = np.column_stack([
        operator.apply_u(np.eye(m)[:, column])
        for column in range(m)
    ])
    laplacian = operator.laplacian.toarray()
    p = u @ u.T
    q = np.eye(n) - p
    b = q @ laplacian @ u
    app = u.T @ laplacian @ u
    h = q @ laplacian @ q + p
    inverse = np.linalg.inv(h)
    correction = b.T @ inverse @ b
    slope = b.T @ inverse @ inverse @ b
    return u, p, q, b, app, correction, slope


def test_fibre_isometry_and_orthogonal_projectors(cycle_operator):
    op = cycle_operator
    u, p, q, *_ = _dense_reference(op)
    assert np.allclose(op.apply_ut(op.apply_u(np.eye(op.n_fibres)[:, 0])), [1.0, 0.0])
    assert np.allclose(op.apply_ut(op.apply_u(np.eye(op.n_fibres)[:, 1])), [0.0, 1.0])
    rng = np.random.default_rng(20261119)
    vector = rng.normal(size=op.n_vertices)
    assert np.allclose(op.apply_p(vector), p @ vector, atol=1e-12)
    assert np.allclose(op.apply_q(vector), q @ vector, atol=1e-12)
    assert np.allclose(p @ p, p, atol=1e-12)
    assert np.allclose(q @ q, q, atol=1e-12)
    assert np.allclose(p.T, p, atol=1e-12)
    assert np.allclose(q.T, q, atol=1e-12)
    assert np.allclose(u.T @ u, np.eye(op.n_fibres), atol=1e-12)


def test_matrix_free_b_and_bt_are_adjoint(cycle_operator):
    op = cycle_operator
    rng = np.random.default_rng(12345)
    coarse = rng.normal(size=op.n_fibres)
    fine = rng.normal(size=op.n_vertices)
    lhs = float(np.dot(op.apply_b(coarse), fine))
    rhs = float(np.dot(coarse, op.apply_bt(fine)))
    assert abs(lhs - rhs) < 1e-12


def test_matrix_free_galerkin_block_matches_dense_reference(cycle_operator):
    op = cycle_operator
    _, _, _, _, app, *_ = _dense_reference(op)
    for coarse in (np.array([1.0, 0.0]), np.array([0.25, -2.0])):
        assert np.allclose(op.apply_app(coarse), app @ coarse, atol=1e-12)


def test_matrix_free_feshbach_correction_matches_dense_schur_complement(cycle_operator):
    op = cycle_operator
    *_, correction, slope = _dense_reference(op)
    rng = np.random.default_rng(2026)
    coarse = rng.normal(size=op.n_fibres)
    actual, info = op.apply_feshbach_correction(coarse)
    assert np.allclose(actual, correction @ coarse, atol=1e-10)
    assert info["info"] == 0
    assert info["residual_norm"] < 1e-10
    assert np.linalg.eigvalsh(correction).min() >= -1e-12


def test_matrix_free_self_energy_slope_matches_dense_second_resolvent(cycle_operator):
    op = cycle_operator
    *_, slope = _dense_reference(op)
    coarse = np.array([0.75, -1.25])
    actual, info = op.apply_self_energy_slope(coarse)
    assert np.allclose(actual, slope @ coarse, atol=1e-10)
    assert info["first_solve"]["residual_norm"] < 1e-10
    assert info["second_solve"]["residual_norm"] < 1e-10
    assert np.linalg.eigvalsh(slope).min() >= -1e-12


def test_matrix_free_operator_rejects_invalid_graph_and_fibres():
    with pytest.raises(ValueError, match="self-loops"):
        MatrixFreeFeshbach(np.array([[0, 0]]), np.array([0]))
    with pytest.raises(ValueError, match="duplicate"):
        MatrixFreeFeshbach(np.array([[0, 1], [1, 0]]), np.array([0, 1]))
    with pytest.raises(ValueError, match="nonempty"):
        MatrixFreeFeshbach(np.array([[0, 1]]), np.array([], dtype=int))
    with pytest.raises(ValueError, match="nonempty"):
        MatrixFreeFeshbach(np.array([[0, 1]]), np.array([0, 2]))
